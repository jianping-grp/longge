"""
Backend for running residue mutation calculations. Used by the
Residue Scanning, Antibody Humanization, and Affinity Maturation
panels.

Properties are added to structure objects as CT-level properties. These are
added as before and after values. Properties available include:
   - potential energy
   - prime energy
   - pKa
   - SASA (polar)
   - SASA (non-polar)
   - SASA (total)
   - hydrophobicty
   - hydrophilicity
   - solubility
   - rotatable bonds
   - vdw surface complementarity

Copyright (c) Schrodinger, LLC. All rights reserved
"""
from __future__ import print_function

#- Imports ------------------------------------------------------------------

from __future__ import absolute_import
from __future__ import division
from past.utils import old_div
import sys
import os
import re
import time
import traceback
import argparse
import errno
import shutil

from schrodinger.application.bioluminate import protein
from schrodinger.job import jobcontrol, launcher, queue
from schrodinger import structure
from schrodinger.structutils import analyze, minimize, build
from schrodinger.utils import subprocess, fileutils
import schrodinger.utils.log as log
from schrodinger.utils import license as biolicense
from schrodinger.protein import remediate, captermini
import future.utils
from six.moves import range
from six.moves import zip

try:
    from schrodinger.application.prime.packages import rottemp
    import schrodinger.application.prime.packages.utilities as psp_util
    from schrodinger.application.prime.packages.ncaa import prep_ncaa
except ImportError:
    rottemp = None
    psp_util = None
    prep_ncaa = None

#- Globals ------------------------------------------------------------------

# these regex support arbitrary 3 letter code for residue names
# by default Mutator class only supports standard residues
MULTI_MUTATION_RE = re.compile(r"""
    (?P<chain>[a-zA-Z_]{1})
    :
    (?P<resnum>-?\d+)
    (?P<inscode>[a-zA-Z]{1})?  # optional
    \s?                        # optional
    (?P<mutations>([a-zA-Z0-9]{3})(,[a-zA-Z0-9]{3})*)? # optional
    """, re.VERBOSE)

# here new_resname can be three letter new residue or loop insertion/deletion
SINGLE_MUTATION_RE = re.compile(r"""
    (?P<chain>[a-zA-Z_]{1})
    :
    (?P<resnum>-?\d+)
    (?P<inscode>[a-zA-Z]{1})?  # optional
    ->
    (?P<new_resname>(\+?\-?[a-zA-Z0-9]+))
    """, re.VERBOSE)

# This regex handles capability to match Mutation strings returned by Prime
# (stored as ct properties)
PRIME_MUTANT_RE = re.compile(r"""
    (?P<chain>[a-zA-z_]{1})
    :
    (?P<resnum>\d+)
    (?P<inscode>[a-zA-Z]{1})?  # optional
    \(
    (?P<old_resname>[a-zA-Z]{3})
    ->
    (?P<new_resname>(\+?\-?[a-zA-Z0-9]+))
    \)
    """, re.VERBOSE)

logger = log.get_output_logger('residue_scanning_backend.py')
logger.setLevel(log.logging.INFO)
if log.get_environ_log_level() <= log.DEBUG:
    logger.setLevel(log.DEBUG)

PRIME = os.path.join(os.environ["SCHRODINGER"], "prime")
PRIMEEXE = os.path.join(os.environ["SCHRODINGER"], "prime.exe")

# All of the properties used in this script
DELTA_PROP_BASE = 'r_bioluminate_delta_'
SUBUNIT_PROP = 's_bioluminate_subunit_asl'
DELTA_POTSTABILITY_PROP = DELTA_PROP_BASE + 'Potential_Stability'  # do we really need this?
DELTA_AFFINITY_PROP = DELTA_PROP_BASE + 'Affinity'
DELTA_STABILITY_PROP = DELTA_PROP_BASE + 'Stability'
FEP_SUITABILITY_PROP = DELTA_PROP_BASE + 'FEP_Suitability'
UNFOLDED_CONTRIBUTION_PROP = protein.Mutator.UNFOLDED_PROPERTY
UNFOLDED_CONTRIBUTION_PRIME_PROP = protein.Mutator.UNFOLDED_PROPERTY_PRIME
TYPE_PROP = 's_bioluminate_Entry_type'

# The types for each protein and its subunits
COMPLEX_TYPE = 'Complex'
RECEPTOR_TYPE = 'Receptor'
LIGAND_TYPE = 'Ligand'

BACKEND = 'residue_scanning_backend.py'

#- Classes ------------------------------------------------------------------


class RefinementError(Exception):
    """ A custom exception for any refinement failures """


class Scanner:
    """
    Class that scans proteins for subunits and handles all mutations and
    property calculations.

    """
    LIGAND_SUBUNIT = 'ligand'
    PROTEIN_SUBUNIT = 'protein'
    NO_SUBUNIT = None

    def __init__(self,
                 reference_st,
                 refine_mut,
                 subunit_type=None,
                 asl=None,
                 solvent='vsgb2.0',
                 use_membrane=False,
                 loop=None,
                 loop_opt=None,
                 loop2=None,
                 cleanup=True,
                 nbcutoff=14.0):
        """
        :param reference_st: The structure used as the reference.
        :param refine_mut: The refinement to carry out on the mutated st.
        :param subunit_type: Type of subunit. See class props above.
        :param asl: The asl used to determine the subunit.
        :param solvent: The solvent model for Prime job.
        :param loop: The loop residues for Prime loop refinement.
        :param loop_opt: The options for Prime loop refinement.
        :param loop2: The second loop residues for Prime loop refinement.
        :param nbcutoff: The non-bonded cutoff for Python minimizer.
        """
        self.original_st = reference_st
        self.refine_mut = refine_mut

        if self.refine_mut.startswith('prime_'):
            self.calculate_prime_energy = True
        else:
            self.calculate_prime_energy = False

        self.subunit_type = subunit_type
        self.subunit_asl = asl

        self.solvent = solvent
        self.use_membrane = use_membrane
        self.loop = loop
        self.loop_opt = loop_opt
        self.loop2 = loop2
        self.cleanup = cleanup
        self._min_nbcutoff = nbcutoff

        self.ref_subunits = self.getSubunits(self.original_st, True)

    def _clearSubunitProperties(self, st):
        """
        Clears the properties from a complex to prepare it to be a subunit.
        these properties need to be removed otherwise incorrect data can
        sneak in.

        """
        st2 = st.copy()
        mutations = st.property.get(protein.Mutator.MUTATIONS_PROPERTY, '')
        st.property.clear()
        st.property[protein.Mutator.MUTATIONS_PROPERTY] = mutations

        # keep implicit membrane
        for prop in st2.property:
            if 'r_psp_Mem' in prop:
                st.property[prop] = st2.property[prop]

    def getSubunits(self, st, is_reference=False):
        """
        Get the references subunits based on the asl supplied in the
        initialization.

        :param is_reference: Whether the subunits we are extracting are
                from the reference st or not.

        """
        # We need to copy the complex sent in and clear all of its properties.
        # These properties can be passed on and give incorrect information.
        # Fix for Ev:133073
        self._clearSubunitProperties(st)

        if not self.subunit_asl:
            return [st]

        # "subunit1" is the receptor, "subunit2" is the ligand
        if self.subunit_type == self.LIGAND_SUBUNIT:
            asl_1 = '(NOT (%s))' % self.subunit_asl
            asl_2 = self.subunit_asl
        elif self.subunit_type == self.PROTEIN_SUBUNIT:
            asl_1 = self.subunit_asl
            asl_2 = '(NOT (%s))' % self.subunit_asl

        atom_list = analyze.evaluate_asl(st, asl_1)
        # Extract the atoms and make sure to copy the properties. There are
        # properties added by the Mutator class that we need.
        subunit1 = st.extract(atom_list, copy_props=True)
        subunit1.property[SUBUNIT_PROP] = asl_1

        atom_list = analyze.evaluate_asl(st, asl_2)
        subunit2 = st.extract(atom_list, copy_props=True)
        subunit2.property[SUBUNIT_PROP] = asl_2

        if self.subunit_type:
            # Add the entry type to the subunits
            subunit1.property[TYPE_PROP] = RECEPTOR_TYPE
            subunit2.property[TYPE_PROP] = LIGAND_TYPE

            if is_reference:
                subunit1.title = '%s (Reference, unrefined receptor)' % st.title
                subunit2.title = '%s (Reference, unrefined ligand)' % st.title
            else:
                subunit1.title = '%s (receptor)' % st.title
                subunit2.title = '%s (ligand)' % st.title

        return [subunit1, subunit2]

    def _refineResidues(self, st, residues, jobname):
        """
        Run the refinement. If the refinement type uses Prime we will get
        many 'r_psp_Prime_<prop>' properties.

        """
        if not residues:
            raise ValueError("Please specify a list of residues to refine")
            # Otherwise Prime will refine the whole CT

        try:
            refiner = protein.Refiner(st, residues=residues)
            if self.use_membrane:
                mem_str = 'yes'
            else:
                mem_str = 'no'
            if self.loop:
                res_sphere = 0.0
                maxcalpha = None
                protocol = 'LOOP_BLD'
                start_res = self.loop[0]
                end_res = self.loop[1]
                loop2 = None
                host = None
                if self.loop_opt:
                    for opt in self.loop_opt:
                        k, v = opt.split('/')
                        if k == 'RES_SPHERE':
                            res_sphere = float(v)
                        elif k == 'MAX_CA_MOVEMENT':
                            maxcalpha = float(v)
                        elif k == 'PROTOCOL':
                            protocol = v
                        elif k == 'HOST':
                            host = v
                        else:
                            raise ValueError(
                                'Non-supported keyword/value: %s' % opt)
                if self.loop2:
                    protocol = 'LOOP_PAIR'
                    loop2 = self.loop2[:2]
                    # if two loops are in differnet subunits, then LOOP_PAIR
                    # does not work
                    if (not residue_in_struct(
                            st, self.loop2[0])) and residue_in_struct(
                                st, self.loop[0]):
                        protocol = 'LOOP_BLD'
                        loop2 = None
                    elif (not residue_in_struct(
                            st, self.loop[0])) and residue_in_struct(
                                st, self.loop2[0]):
                        protocol = 'LOOP_BLD'
                        loop2 = None
                        start_res = self.loop2[0]
                        end_res = self.loop2[1]
                    elif (not residue_in_struct(st, self.loop[0])) and (
                            not residue_in_struct(st, self.loop2[0])):
                        raise ValueError(
                            'either loop or loop2 is in the structure!')
                refined_st = refiner.runRefinement(
                    self.refine_mut,
                    jobname,
                    SGB_MOD=self.solvent,
                    USE_MEMBRANE=mem_str,
                    start_res=start_res,
                    end_res=end_res,
                    res_sphere=res_sphere,
                    maxcalpha=maxcalpha,
                    protocol=protocol,
                    loop2=loop2,
                    host=host)
            else:
                refined_st = refiner.runRefinement(
                    self.refine_mut,
                    jobname,
                    SGB_MOD=self.solvent,
                    USE_MEMBRANE=mem_str)
            # All arguments except the first two are passed to Prime

            if self.cleanup:
                refiner.clean()
        except Exception as err:
            msg = '\nException: %s' % str(err)
            msg += '\nSkipping this mutation. Problems encountered during refinement. '
            msg += 'Please make sure the protein is properly prepared.'
            raise RefinementError(msg)

        # BIOLUM-1315 Let the Viewer GUI know which refinement method we used:
        refined_st.property["s_bioluminate_Refinement_method"] = self.refine_mut

        return refined_st

    def _refineAndCalculateProps(self, st, calculations, jobname, refine=True):
        """
        Get the properties for a structure

        """
        # Prime is picky and if you ask it to refine only a subset of residues
        # each in the subset needs to be there or it will fail
        refine_residues = [
            r for r in st.residue if str(r) in self.refine_residues
        ]

        residues_str = ", ".join([str(res) for res in refine_residues])
        info("  Refining residues: %s" % residues_str)

        # Perform refinement only when requested. Note that if refinement is not
        # performed, then the structure will not have a "Prime Energy" property.
        # See BIOLUM-1423 (a fix for which was reverted).
        t1 = time.time()
        if refine:
            refine_res_string = ", ".join([str(r) for r in refine_residues])
            debug("Refining residues: %s" % refine_res_string)

            refine_jobname = 'refine_job_%s' % jobname
            st = self._refineResidues(st, refine_residues, refine_jobname)
        t2 = time.time()
        info('TIME for refinement %f seconds' % (t2 - t1))

        props = calculate_props(
            st,
            jobname,
            calculations,
            cleanup=self.cleanup,
            nbcutoff=self._min_nbcutoff,
            calculate_prime_energy=self.calculate_prime_energy,
            lig_asl=self.subunit_asl,
        )
        t3 = time.time()
        info('TIME for prop calc %f seconds' % (t3 - t2))
        return (st, props)

    def _setFinalDeltaEnergyProp(self, st, include_stability):
        # NOTE: The unfolded contribution is already a "delta":
        unfolded_contribution_gas = st.property.get(UNFOLDED_CONTRIBUTION_PROP,
                                                    0.0)
        unfolded_contribution_prime = st.property.get(
            UNFOLDED_CONTRIBUTION_PRIME_PROP, 0.0)

        # Get the energy property based on refinement type. After the
        # subunits are calculated we change it to "r_bioluminate_delta_Potential_Stability"
        # and we add "r_bioluminate_delta_Affinity" to the complex for affinity
        # calculations and "r_bioluminate_delta_Stability" for stability
        # calculations.

        # NOTE: This is why we are forcing e_pot to always be calculated in
        # -classic mode:
        e_pot_prop = DELTA_PROP_BASE + protein.PROP_CONVERSIONS['e_pot']
        st.property[
            DELTA_POTSTABILITY_PROP] = st.property[e_pot_prop] - unfolded_contribution_gas
        # NOTE: This is correct. The delta energy is the Minimizer energy minus unfolded contribution
        # See BIOLUM-1267

        if include_stability:
            if self.calculate_prime_energy:
                # Prime was used for refinement:
                prime_e_prop = DELTA_PROP_BASE + protein.PROP_CONVERSIONS['prime_energy']
                st.property[
                    DELTA_STABILITY_PROP] = st.property[prime_e_prop] - unfolded_contribution_prime
            else:
                try:
                    # The python minimization was used for refinement:
                    st.property[
                        DELTA_STABILITY_PROP] = st.property[e_pot_prop] - unfolded_contribution_gas
                    # NOTE: delta stability is equal to delta potential energy
                    # in this case
                except:
                    # otherwise, set to 0.0
                    st.property[DELTA_STABILITY_PROP]

    def _getEnergy(self, st):
        """
        Returns the energy for this structure. If Prime was used for refinement,
        returns the value of the "prime_energy" property; otherwise returns the
        value of the "e_pot" property.
        """
        if self.calculate_prime_energy:
            energy_prop = DELTA_PROP_BASE + protein.PROP_CONVERSIONS['prime_energy']
        else:
            energy_prop = DELTA_PROP_BASE + protein.PROP_CONVERSIONS['e_pot']
        return st.property[energy_prop]

    def _calculateAffinity(self, calculations, jobname):
        """
        Run the supplied `calculations` for:
           - The reference protein
           - The mutated protein
           - Subunit 1 of the reference protein
           - Subunit 1 of the mutated protein
           - Subunit 2 of the reference protein
           - Subunit 2 of the mutated protein

        """
        info("\nCalculating properties for refined reference structure...")
        refined_ref, ref_agg_props = self._refineAndCalculateProps(
            self.reference.copy(), calculations, 'ref_' + jobname)

        info("Calculating properties for refined mutant structure...")
        mut_st, mut_agg_props = self._refineAndCalculateProps(
            self.mutant.copy(), calculations, 'mut_' + jobname)

        calculated_sts = []

        # Grab the aggregate calculations and apply them to the mutated
        # structure
        set_delta_props(mut_st, ref_agg_props, mut_agg_props)
        mut_st.property[TYPE_PROP] = COMPLEX_TYPE

        calculated_sts.append(mut_st)

        ref_subunits = self.getSubunits(self.reference)
        mut_subunits = self.getSubunits(self.mutant)

        # Grab the subunit calculations and apply them to the mutated subunit
        # structures
        sub_energies = 0.0

        # Change the property names to generic Energy property for
        # python- and prime-specific energy

        self._setFinalDeltaEnergyProp(mut_st, include_stability=True)

        # The binding affinity will have its subunit's energies subtracted
        # OLD: affinity = mut_st.property[DELTA_POTSTABILITY_PROP]
        affinity = self._getEnergy(mut_st)
        debug("  Delta energy of the full structure: %.2f" %
              self._getEnergy(mut_st))

        for ref, mut in zip(ref_subunits, mut_subunits):
            stype = ref.property[TYPE_PROP]  # be either ligand or receptor
            info("Calculating properties for subunit reference %s structure..."
                 % stype)

            # Test just in case:
            if ref.property[TYPE_PROP] != mut.property[TYPE_PROP]:
                raise ValueError("reference and mutant subunit type mismatch")

            # Figure out where the mutation residues are, on ligand or receptor
            # and skip the refinement if the subunit does not have mutation
            # residues
            mutate_residues = [
                r for r in ref.residue if str(r) in self.mutate_residues
            ]
            if mutate_residues:
                refine = True
            else:
                refine = False

            ref_sub_st, ref_sub_props = self._refineAndCalculateProps(
                ref, calculations, 'ref_sub_' + jobname, refine)

            stype = mut.property[TYPE_PROP]
            info("Calculating properties for subunit mutated %s structure..." %
                 stype)

            mut_sub_st, mut_sub_props = self._refineAndCalculateProps(
                mut, calculations, 'mut_sub_' + jobname, refine)

            set_delta_props(mut_sub_st, ref_sub_props, mut_sub_props)

            # Change the property names to generic Energy property for
            # python- and prime-specific energy

            self._setFinalDeltaEnergyProp(mut_sub_st, include_stability=False)

            calculated_sts.append(mut_sub_st)

            # Subtract the subunit energies from affinity
            # OLD: affinity -= mut_sub_st.property[DELTA_POTSTABILITY_PROP]
            affinity -= self._getEnergy(mut_sub_st)
            debug("  Delta energy of the %s subunit structure: %.2f" %
                  (mut.property[TYPE_PROP], self._getEnergy(mut_sub_st)))

        mut_st.property[DELTA_AFFINITY_PROP] = affinity
        debug("  Calculated delta affinity: %.2f" % affinity)

        return refined_ref, calculated_sts

    #       ^
    # FIXME    Combine these methods?  Or at factor out shared code?
    #       v

    def _calculateStability(self, calculations, jobname):
        """
        Run the supplied `calculations` for a reference protein and its
        mutant.

        """
        info("\nCalculating properties for the refined reference structure...")
        refined_ref, ref_agg_props = self._refineAndCalculateProps(
            self.reference.copy(), calculations, 'ref_' + jobname)

        info("Calculating properties for the refined mutant structure...")
        mut_st, mut_agg_props = self._refineAndCalculateProps(
            self.mutant.copy(), calculations, 'mut_' + jobname)

        set_delta_props(mut_st, ref_agg_props, mut_agg_props)

        self._setFinalDeltaEnergyProp(mut_st, include_stability=True)

        mut_st.property[TYPE_PROP] = COMPLEX_TYPE

        return refined_ref, [mut_st]

    def _setRefineResidues(self, residue_map, dist):
        """
        Sets the `self.refine_residues` which is used to later calculate
        the properties.
        """

        # Make a list of residue names that will be mutated
        residues_to_mutate = [
            str(r) for r in future.utils.listvalues(residue_map)
        ]

        self.mutate_residues = residues_to_mutate
        self.refine_residues = self.getResiduesWithinMutantArg(
            self.reference, residues_to_mutate, dist)

    def getResiduesWithinMutantArg(self, struct, residues, distance):
        """
        In order to get a standard set of residues to refine for all mutations,
        for each mutation in question, we mutate (again) the residue to an Arg
        and then get the residues within distance of the Arg.  This helps make
        sure that the same set of residues is always refined.

        :type struct: `schrodinger.structure.Structure`
        :param struct: The structure to analyze

        :type residues: list
        :param residues: A list of mutatable residue names  (e.g. "A:28")

        :type distance: float
        :param distance: Find residues within this distance of the Arg mutation

        :rtype: list
        :return: List of residue names (e.g. "A:29b") to optimize
        """

        ref_atoms_to_refine = set()

        # Set of names of residues to optimize:
        optimize_residues = set()

        for res_str in residues:
            mutant_st = struct.copy()

            # Find an atom from this residue in the given structure:
            ref_atom = mutant_st.findResidue(res_str).getAtomIndices()[0]

            # Mutate this residue to arginine (largets of the 20):
            atom_map = build.mutate(mutant_st, ref_atom, 'ARG')
            # atom_map is a dict that translates old atom numbers to new

            # Find the residue again, since atom numbering has changed:
            new_res = mutant_st.findResidue(res_str)

            residue_asl = protein.get_residue_asl(new_res)
            within_asl = 'within %.2f (%s)' % (distance, residue_asl)
            # atoms_within_mut are numbered via the mutated structure numbering
            atoms_within_mut = set(analyze.evaluate_asl(mutant_st, within_asl))

            for anum in atoms_within_mut:
                res = mutant_st.atom[anum].getResidue()
                optimize_residues.add(str(res))

        return optimize_residues

    def calculateOneMutation(self, mutation, dist, calculations, jobname):
        """
        Refine the given mutation CT and calculate properties. The reference structure
        is also self-mutated and refined in exactly the same way for comparison.
        """

        # Set our mutant st and idealized reference st:
        self.mutant = mutation.struct
        self.reference = mutation.ref_struct
        residue_map = mutation.residue_map

        # Set the residues we are going to refine after the calculations are
        # complete
        self._setRefineResidues(residue_map, dist)

        # Run a calculation for protein-protein binding affinity
        # (Gets executed when -receptor_asl argument was specified)
        if self.subunit_type:
            return self._calculateAffinity(calculations, jobname)

        # Run a protein stability calculation
        # (Neither -ligand_asl nor -receptor_asl were specified)
        else:
            return self._calculateStability(calculations, jobname)


class PrimeResidueScan:
    """
    Class that runs a Prime residue mutation job
    """

    def __init__(self,
                 orig_st,
                 mutations_list,
                 jobname,
                 nojobid=False,
                 refine_mut='prime_residue',
                 ligand_asl=None,
                 residue_structure=None,
                 dist=0.0,
                 mutable_only=False,
                 all_mutable=False,
                 rigidmove=False,
                 refine_unbound=False,
                 solvent='vsgb2.0',
                 use_membrane=False,
                 mc_job=False,
                 nstep=2000,
                 noutput=100,
                 temperature=298.0,
                 write_trj=False,
                 nterm=1000000,
                 stab_weight=0.0,
                 aff_weight=1.0,
                 stab_cutoff=30.0,
                 aff_cutoff=30.0,
                 concurrent=1,
                 random_seed=0,
                 random_start=False,
                 cleanup=True,
                 fep_preparation=False):
        """
        :param orig_st: The input structure
        :param mutations_list: A list of mutations ` L{(chain, resnum, inscode, new_resname)` }
        :param dist: The distance cutoff for nearby residues to be refined.
        :param solvent: The solvent model for Prime jobs.
        :param mc_job: If this is a MC residue mutation.
        :param random_seed: The random seed for MC sampling.
        :param random_start: Start the MC search from a random sequence instead of the starting sequence
        :param fep_preparation: Run in FEP preparation mode.
                This returns different output names: (FEP_SUITABILITY_PROP) and
                requires FEP licenses instead of Bioluminate ones (still requires
                PSP_PLOP tokens).
        """

        remediate.remediate_ct(orig_st)  # atom names must be correct
        start_st = orig_st.copy()

        # Prepare for Prime input: 1. annotate mutable residues
        mutable_residues = []
        nonstan = []
        nonstan_cts = []
        if residue_structure:
            for ct in structure.StructureReader(residue_structure):
                nonstan_cts.append(ct.copy())
        for changes in mutations_list:
            for chain, resnum, inscode, new_resname in changes:
                res_id = (chain, resnum, inscode)
                if res_id not in mutable_residues:
                    mutable_residues.append(res_id)
                if new_resname not in rottemp.PlopRotamerTemplateAssignment.standard_residues:
                    if new_resname not in nonstan:
                        nonstan.append(new_resname)
        tail_mut = False  # terminal residues need different templates
        for chain, resnum, inscode in mutable_residues:
            res_str = '%s:%s%s' % (chain, resnum, inscode)
            res = start_st.findResidue(res_str)
            if captermini.res_can_be_capped(res) != 0:
                tail_mut = True
            res_atoms = res.getAtomIndices()
            for iatom in res_atoms:
                start_st.atom[iatom].property['i_psp_Mutable_Residue'] = 1
            # check if the starting structure is a nonstandard residue
            old_resname = start_st.atom[res_atoms[0]].pdbres.replace(' ', '')
            if old_resname not in rottemp.PlopRotamerTemplateAssignment.standard_residues:
                nonstan.append(old_resname)
                res_ct = start_st.extract(res_atoms)
                res_ct.title = old_resname
                nonstan_cts.append(res_ct)
        # Prepare for Prime input: 2. create rot library for non-standard amino
        # acids
        nonstan_data = []
        if nonstan:
            for new_resname in nonstan:
                ncaa_exists = False
                for ct in nonstan_cts:
                    found = False
                    for atom in ct.atom:
                        if atom.pdbres.strip() == new_resname:
                            found = True
                            break
                    if not found:
                        continue
                    ncaa_exists = True
                    template_cts = []
                    rotlib_lines = []
                    # already fully prepared
                    if 's_psp_NCAA_Rotlib' in ct.property and 'r_psp_NCAA_Energy' in ct.property:
                        prepared_ct = ct.copy()
                    else:
                        try:
                            prepared_ct = prep_ncaa(
                                ct,
                                new_resname,
                                position='mid',
                                energy_calc=True,
                                solvent=solvent,
                                cleanup=cleanup)
                        except RuntimeError as err:
                            info(
                                'ERROR: cannot prepare non-standard amino acids %s!'
                                % new_resname)
                            info('\nError message from the preparation:\n' +
                                 str(err))
                            sys.exit(1)
                    rotlib_str = prepared_ct.property['s_psp_NCAA_Rotlib']
                    # possible to not have a library, e.g. ALA
                    if rotlib_str:
                        rotlib_lines = rotlib_str.split(',')
                    else:
                        rotlib_lines = []
                    template_cts.append(prepared_ct)
                    if tail_mut:
                        try:
                            # just need annotated CT, not the rot lib. There
                            # might be a better way.
                            prepared_ct = prep_ncaa(
                                ct,
                                new_resname,
                                position='nterm',
                                energy_calc=False)
                            template_cts.append(prepared_ct)
                            prepared_ct = prep_ncaa(
                                ct,
                                new_resname,
                                position='cterm',
                                energy_calc=False)
                            template_cts.append(prepared_ct)
                        except RuntimeError as err:
                            info(
                                'ERROR: cannot prepare terminal non-standard amino acids %s!'
                                % new_resname)
                            info('\nError message from the preparation:\n' +
                                 str(err))
                            sys.exit(1)
                    # ready to pass to Prime backend
                    nonstan_data.append((new_resname, rotlib_lines,
                                         template_cts))
                if not ncaa_exists:
                    raise RuntimeError(
                        'ERROR: non-standard amino acid is not defined: %s!' %
                        new_resname)

        self.orig_st = orig_st
        self.start_st = start_st
        self.nonstan_data = nonstan_data
        self.jobname = jobname
        self.nojobid = nojobid
        self.refine_mut = refine_mut
        self.dist = dist
        self.mutable_only = mutable_only
        self.all_mutable = all_mutable
        self.rigidmove = rigidmove
        self.refine_unbound = refine_unbound
        self.solvent = solvent
        self.use_membrane = use_membrane
        self.mc_job = mc_job
        self.nstep = nstep
        self.noutput = noutput
        self.temperature = temperature
        self.write_trj = write_trj
        self.nterm = nterm
        self.stab_weight = stab_weight
        self.aff_weight = aff_weight
        self.stab_cutoff = stab_cutoff
        self.aff_cutoff = aff_cutoff
        self.concurrent = concurrent
        self.random_seed = random_seed
        self.random_start = random_start
        self.ligand_asl = ligand_asl
        self.mutations_list = mutations_list
        self.cleanup = cleanup
        self.fep_preparation = fep_preparation

    def _createJobDir(self):
        """
        Create a directory to run residue mutation job
        """
        dirs = [d for d in os.listdir(os.getcwd()) if os.path.isdir(d)]
        count = 0
        for d in dirs:
            match = re.search(r'scan_job_(\d+)', d)
            if match:
                _count = int(match.group(1))
                if _count > count:
                    count = _count
        retries = 5
        while True:
            temp_dir = 'scan_job_%d' % (count + 1)
            try:
                os.mkdir(temp_dir)
                break
            except OSError as err:
                if err.errno != errno.EEXIST:
                    raise
                count += 1
                retries -= 1
                if retries < 1:
                    raise
        return os.path.abspath(temp_dir)

    def clean(self):
        try:
            shutil.rmtree(self.jobdir)
        except Exception as err:
            print('Unable to remove the directory %s\n%s' % (self.jobdir, err))

    def _writeInput(self):
        self.jobdir = self._createJobDir()

        pjobname = self.jobname + '_scan'

        structfile = os.path.join(self.jobdir, pjobname + '_start.maegz')
        self.start_st.write(structfile)
        origfile = os.path.join(self.jobdir, pjobname + '_orig.maegz')
        self.orig_st.write(origfile)

        mutfile = os.path.join(self.jobdir, pjobname + '_mut.txt')
        fh = open(mutfile, 'w')
        for i, mut_changes in enumerate(self.mutations_list, start=1):
            mut_strings = []
            for chain, resnum, inscode, new_resname in mut_changes:
                if inscode == " ":
                    inscode = ""
                mut_str = "%s:%s%s %s" % (chain, resnum, inscode, new_resname)
                mut_strings.append(mut_str)
            fh.write(" ".join(mut_strings) + "\n")
        fh.close()

        inpfile = os.path.join(self.jobdir, pjobname + '.inp')
        fhinp = open(inpfile, 'w')
        fhinp.write("""\
STRUCT_FILE     %s
NATIVE_FILE     %s
MUTATION_FILE   %s
PRIME_TYPE      RESIDUE_SCAN
OPLS_VERSION    2005
""" % (structfile, origfile, mutfile))

        matom_nonstan = 0
        if self.nonstan_data:
            residue_sfile = os.path.join(self.jobdir,
                                         pjobname + '_nonstan_struct.mae')
            if os.path.exists(residue_sfile):
                os.remove(residue_sfile)
            for i, (resname, lib_contents, cts) in enumerate(self.nonstan_data):
                for ct in cts:
                    ct.append(residue_sfile)
                    if len(ct.atom) > matom_nonstan:
                        matom_nonstan = len(ct.atom)
                # if non-empty, write a side lib file for PLOP to use
                if lib_contents:
                    residue_rotfile = os.path.join(self.jobdir,
                                                   resname + '_01.side')
                    fh = open(residue_rotfile, 'w')
                    for line in lib_contents:
                        fh.write(line)
                    fh.close()
            fhinp.write('NONSTAN_RESIDUE_FILE %s\n' % residue_sfile)
            fhinp.write('NONSTAN_ROT_DIR %s\n' % self.jobdir)
        matom_mutres = max(matom_nonstan, 24) + 20  # add enough buffer
        fhinp.write('MAX_MUTABLE_RESIDUE_SIZE %d\n' % matom_mutres)

        if self.refine_mut.startswith('prime_'):
            fhinp.write('SCAN_REFINE_METHOD %s\n' % self.refine_mut[6:])
        if self.dist > 0.0:
            fhinp.write('SCAN_DISTANCE_CUTOFF %f\n' % self.dist)
        if self.mutable_only:
            fhinp.write('SCAN_MUTABLE_ONLY yes\n')
        if self.all_mutable:
            fhinp.write('SCAN_ALL_MUTABLE yes\n')
        if self.rigidmove:
            fhinp.write('SCAN_RIGIDBODY_MOVE yes\n')
        if self.refine_unbound:
            fhinp.write('SCAN_REFINE_UNBOUND yes\n')
        # only write out non default energy model
        if self.solvent != 'vsgb2.0':
            fhinp.write('SGB_MOD  %s\n' % self.solvent)
        if self.use_membrane:
            fhinp.write('USE_MEMBRANE yes \n')
        if self.ligand_asl:
            fhinp.write("LIGAND          asl = (%s) \n" % self.ligand_asl)
        if self.mc_job:
            fhinp.write("SCAN_TYPE MC\n")
            fhinp.write("SCAN_CONCURRENT  %s \n" % self.concurrent)
            fhinp.write("SCAN_NSTEP  %i \n" % self.nstep)
            fhinp.write("SCAN_NOUTPUT  %i \n" % self.noutput)
            fhinp.write("SCAN_TEMP  %f \n" % self.temperature)
            if self.write_trj:
                fhinp.write("SCAN_WRITE_TRJ  yes\n")
            fhinp.write("SCAN_NTERM %i \n" % self.nterm)
            fhinp.write("SCAN_WSTAB %f \n" % self.stab_weight)
            fhinp.write("SCAN_WAFF %f \n" % self.aff_weight)
            fhinp.write("SCAN_STAB_CUTOFF %f \n" % self.stab_cutoff)
            fhinp.write("SCAN_AFF_CUTOFF %f \n" % self.aff_cutoff)
            # does not support arbitrary combination yet
            if self.stab_weight >= 0.99 and self.aff_weight <= 0.01:
                fhinp.write("OPT_PROPERTY    stability \n")
            elif self.aff_weight >= 0.99 and self.stab_weight <= 0.01:
                fhinp.write("OPT_PROPERTY    affinity \n")
            fhinp.write("SCAN_RANDOM_SEED  %i \n" % self.random_seed)
            if self.random_start:
                fhinp.write("SCAN_RANDOM_START  yes\n")
        fhinp.close()
        self.inpfile = inpfile

    def run(self):
        self._writeInput()
        init_dir = os.getcwd()
        os.chdir(self.jobdir)

        command = ["prime", os.path.basename(self.inpfile)]
        if self.nojobid:
            command.append('-NOJOBID')

        # Run the command
        if not self.nojobid:
            try:
                job = jobcontrol.launch_job(command)
            except Exception as msg:
                os.chdir(init_dir)
                raise RuntimeError('Job failed to launch.\nJob file: %s'
                                   '\n\nError Returned:\n%s\n' % (self.inpfile,
                                                                  str(msg)))
            os.chdir(init_dir)
            job.wait()

            if not job.succeeded():
                lastlines = []
                logfile = fileutils.strip_extension(
                    self.inpfile) + '.log'  # replace .inp with .log
                with open(logfile) as fh:
                    alines = fh.readlines()
                    lastlines = alines[-30:]
                errmsg = ''.join(lastlines)
                raise RuntimeError(
                    'Prime Residue Scanning Job failed.\nJob file: %s\nJobId: %s\n'
                    '\nBelow are the last 30 lines of the log file...\n\n%s' %
                    (self.inpfile, job.job_id, errmsg))
        else:  # nojobid
            try:
                outp = subprocess.check_output(command, universal_newlines=True)
            except subprocess.CalledProcessError as err:
                lastlines = outp.split('\n')[-30:]
                raise RuntimeError(
                    'Prime Residue Scanning calculation failed.\n'
                    'with error "%s".\n'
                    'Below are the last 30 lines of the output...\n\n%s' %
                    (err, lastlines))
            os.chdir(init_dir)
        # set expected output file names
        # Dev note: 7/17/17
        # this is required since PXJob truncates the job at 60 characters
        # PXJobnames truncated around here:
        # http://opengrok.schrodinger.com/xref/psp-src/python/modules/PXjobtypes.py#244
        # mutants and mutrefs filenames set around here:
        # http://opengrok.schrodinger.com/xref/psp-src/python/modules/PXinput.py#2834
        outjobfiles = os.listdir(self.jobdir)
        self.outfiles = [[
            os.path.join(self.jobdir, f)
            for f in outjobfiles
            if f.endswith('.mutants.maegz')
        ][0]]
        self.outfiles.append([
            os.path.join(self.jobdir, f)
            for f in outjobfiles
            if f.endswith('.mutrefs.maegz')
        ][0])

    def calcDeltaProp(self, calculations, min_nbcutoff, refs_file, out_file):
        info('\nReferences (refined structures) is being written to: %s' %
             refs_file)
        ref_props_all = {}
        for st in structure.StructureReader(self.outfiles[1]):
            st.append(refs_file)
            try:
                desc = st.property['s_psp_Prime_Mutations']
            except KeyError:
                continue

            # skip the starting structure without mutations
            if desc == 'NONE':
                continue

            try:
                props = calculate_props(
                    st,
                    self.jobname,
                    calculations,
                    cleanup=self.cleanup,
                    nbcutoff=min_nbcutoff,
                    calculate_prime_energy=True,
                    lig_asl=self.ligand_asl,
                    fep_preparation=self.fep_preparation,
                )
            except:  # FIXME: Is there a known exception we're trying to catch?
                traceback.print_exc()
                props = []
                continue

            if (desc == 'MINIMIZED'):  # minimized starting structure for MC job
                ref_props_all[desc] = props
            else:
                residues = ''
                for onemut in desc.split(','):
                    residues += onemut[0:(onemut.find('(') + 1)]
                ref_props_all[residues] = props

        info('\nOutput is being written to: %s' % out_file)
        num_out = 0
        for st in structure.StructureReader(self.outfiles[0]):
            try:
                mut_props = calculate_props(
                    st,
                    self.jobname,
                    calculations,
                    cleanup=self.cleanup,
                    nbcutoff=min_nbcutoff,
                    calculate_prime_energy=True,
                    lig_asl=self.ligand_asl,
                    fep_preparation=self.fep_preparation,
                )
            except Exception as err:
                warning("WARNING: Failed to calculate properties on mutant: %s"
                        % st.property['s_psp_Prime_Mutations'])
                warning("Exception message: %s" % err.message)
                mut_props = []

            # need to change the mutation property name from Prime to
            # Bioluminate
            mutkey = 's_psp_Prime_Mutations'
            desc = st.property[mutkey]
            st.property['s_bioluminate_Mutations'] = desc
            del st.property[mutkey]

            # skip the un-mutated structures
            if desc == 'NONE' or desc == 'MINIMIZED':
                st.append(out_file)
                continue

            if mut_props:
                residues = ''
                for onemut in desc.split(','):
                    residues += onemut[0:(onemut.find('(') + 1)]
                if residues in ref_props_all:
                    ref_props = ref_props_all[residues]
                else:
                    try:
                        ref_props = ref_props_all['MINIMIZED']
                    except KeyError:
                        error(
                            'ERROR: cannot find the reference state for mutation: '
                            + desc)
                        continue
                # get the "deltas"
                set_delta_props(st, ref_props, mut_props)

            # now set delta pot-stability, stability, affinity
            # 1. calculate Potential stability
            if 'e_pot' in calculations:
                unfolded_contribution_gas = 0.0
                for onemut in desc.split(','):
                    index1 = onemut.find(
                        '(')  # the format is like "A:5(ASP->GLN)"
                    index2 = onemut.find('>')
                    index3 = onemut.find(')')
                    orig_resname = onemut[index1 + 1:index2 - 1]
                    new_resname = onemut[index2 + 1:index3]
                    try:
                        dG = protein.Mutator.GXG_DATA[new_resname] - protein.Mutator.GXG_DATA[orig_resname]
                    except:
                        dG = 0.0
                    unfolded_contribution_gas += dG
                st.property[
                    UNFOLDED_CONTRIBUTION_PROP] = unfolded_contribution_gas
                e_pot_prop = DELTA_PROP_BASE + protein.PROP_CONVERSIONS['e_pot']
                st.property[
                    DELTA_POTSTABILITY_PROP] = st.property[e_pot_prop] - unfolded_contribution_gas

            # 2. copy over Prime stability
            if self.fep_preparation:
                st.property[FEP_SUITABILITY_PROP] = \
                    st.property['r_psp_Prime_delta_Stability']
            else:
                st.property[DELTA_STABILITY_PROP] = \
                    st.property['r_psp_Prime_delta_Stability']

                # 3. affinity may or may not exist
                try:
                    st.property[DELTA_AFFINITY_PROP] = \
                        st.property['r_psp_Prime_delta_Affinity']
                except KeyError:
                    pass

            num_out += 1
            if not st.title:
                st.title = 'Mutant_' + str(num_out)
            else:
                st.title = st.title + '_Mutant_' + str(num_out)
            st.append(out_file)
        info('\nNumber of mutants generated: %i' % num_out)


#- Functions ----------------------------------------------------------------


def residue_in_struct(ct, resid):  # should move this to psp_util later
    # resid in the format like A:10B or C:9
    found = False
    for res in ct.residue:
        if resid == str(res):
            found = True
            break
    return found


def error(text):
    logger.error(text)


def warning(text):
    logger.warning(text)


def info(text):
    logger.info(text)


def debug(text):
    logger.debug(text)


def calc_items_per_job(total_items, njobs):
    """
    Determine the number of items that should go into each subjob
    """

    items_per_job = int(old_div(float(total_items), float(njobs)) + 0.6)

    # Prevent creation of extra tiny subjob at end:
    while items_per_job * njobs < total_items:
        items_per_job += 1

    return items_per_job


def chunks(l, n):
    """
    Yield successive n-sized chunks from list <l>.
    """
    for i in range(0, len(l), n):
        yield l[i:i + n]


def block_size(t, n):
    """
    Return block sizes of range 1-t divided by n
    The last block may be smaller than others
    """
    start_indices = [i for i in range(0, t, n)]
    for i, index in enumerate(start_indices):
        if i < len(start_indices) - 1:
            yield (start_indices[i + 1] - start_indices[i])
        else:
            yield (t - start_indices[i])


def set_delta_props(st, ref_props, mut_props):
    """
    Processes all the changes from the reference to the mutant. For each
    property calculated a delta value will be added to the structure. The
    property will be <DELTA_PROP_BASE>_<property> where property is a valid
    calculation type from the calc type in the parser.
    """
    for k, ref_val in future.utils.listitems(ref_props):
        prop = DELTA_PROP_BASE + protein.PROP_CONVERSIONS[k]
        mut_val = mut_props.get(k)
        if mut_val is not None:
            st.property[prop] = mut_val - ref_val


def calculate_props(
        st,
        jobname,
        calculations,
        cleanup=True,
        nbcutoff=14.0,
        calculate_prime_energy=False,
        lig_asl=None,
        fep_preparation=False,
):
    """
    A simple function to calculate all requested properties and also handle
    Prime related issues.
    """
    calc = protein.PropertyCalculator(
        st, jobname, cleanup=cleanup, nbcutoff=nbcutoff, lig_asl=lig_asl)
    # Catch propka failures and ignore them, rerun without propka.
    # Print the error to stdout so the log picks it up
    try:
        props = calc.calculate(*calculations)
    except protein.PropkaError as err:
        msg = 'ERROR: PROPKA failed, no pKa will be recorded for this '
        msg += 'calculation.'
        msg += '\n%s' % err
        error(msg)
        calculations = [c for c in calculations if not c == 'pka']
        props = calc.calculate(*calculations)

    if calculate_prime_energy:
        # FIXME this needs re-factoring to remove duplication
        props['prime_energy'] = st.property.get('r_psp_Prime_Energy', 0.0)
        # energy components
        if 'r_psp_Affinity_Covalent' in st.property and not fep_preparation:
            props['affinity_cov'] = st.property.get('r_psp_Affinity_Covalent',
                                                    0.0)
            props['affinity_vdw'] = st.property.get('r_psp_Affinity_vdW', 0.0)
            props['affinity_coulomb'] = st.property.get(
                'r_psp_Affinity_Coulomb', 0.0)
            props['affinity_gb'] = st.property.get('r_psp_Affinity_Solv_GB',
                                                   0.0)
            props['affinity_sa'] = st.property.get('r_psp_Affinity_Solv_SA',
                                                   0.0)
            props['affinity_lipo'] = st.property.get('r_psp_Affinity_Lipo', 0.0)
            props['affinity_hbond'] = st.property.get('r_psp_Affinity_Hbond',
                                                      0.0)
            props['affinity_packing'] = st.property.get(
                'r_psp_Affinity_Packing', 0.0)
            props['affinity_selfcont'] = st.property.get(
                'r_psp_Affinity_SelfCont', 0.0)
        if 'r_psp_Stability_Covalent' in st.property:
            props['stability_cov'] = st.property.get('r_psp_Stability_Covalent',
                                                     0.0)
            props['stability_vdw'] = st.property.get('r_psp_Stability_vdW', 0.0)
            props['stability_coulomb'] = st.property.get(
                'r_psp_Stability_Coulomb', 0.0)
            props['stability_gb'] = st.property.get('r_psp_Stability_Solv_GB',
                                                    0.0)
            props['stability_sa'] = st.property.get('r_psp_Stability_Solv_SA',
                                                    0.0)
            props['stability_lipo'] = st.property.get('r_psp_Stability_Lipo',
                                                      0.0)
            props['stability_hbond'] = st.property.get('r_psp_Stability_Hbond',
                                                       0.0)
            props['stability_packing'] = st.property.get(
                'r_psp_Stability_Packing', 0.0)
            props['stability_selfcont'] = st.property.get(
                'r_psp_Stability_SelfCont', 0.0)
            props['stability_reference'] = st.property.get(
                'r_psp_Stability_Reference', 0.0)

    return props


def get_parser():

    desc = """The driver script to mutate residues and compute various property changes such
    as stability, binding affinity, surface area, etc. It can also search the best mutations
    for affinity maturation with a Monte Carlo based approach. The mutant structures and related
    properties are stored in the output file."""

    refine_choices = (
        protein.Refiner.PYTHON_MINIMIZE,
        protein.Refiner.PRIME_MINIMIZE,  # obsolete, same as PRIME_RESIDUE
        protein.Refiner.PRIME_RESIDUE,
        protein.Refiner.PRIME_SIDECHAIN,
        protein.Refiner.PRIME_SIDECHAIN_CBETA,
        protein.Refiner.PRIME_SIDECHAIN_BB,
        protein.Refiner.PRIME_LOOP_PRED)

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument(
        'structures',
        help='The input structure file containing the receptor, the '
        'receptor/ligand complex, or a Pose Viewer file (first CT is the '
        'receptor, all remaining are docked ligands). If multiple CTs are '
        'in the "struct_file" it is assumed to be a PV file. ')

    parser.add_argument(
        '-jobname',
        default='residue_mutation',
        metavar='',
        help="The base jobname.")
    parser.add_argument(
        '-subjob',
        metavar='',
        default=None,
        help="The number of this subjob (for internal use only)",
    )

    # Specify mutations as source residues and destination types:
    # all residues defined this way have the same destination types
    parser.add_argument(
        '-residues',
        type=residue_type,
        help='List of residues to mutate. Each residue should be in the form of '
        'comma-separated <chain>:<resnum>. For blank chain name, use "_".'
        'Example: "-residues A:122,_:12,A:18A".')
    parser.add_argument(
        '-mut',
        type=mutation_type,
        help='List of residues to mutate reference residues to. Residues should '
        'be a comma-separated list of 3-letter residue names with no spaces.'
        'Example: "-mut ARG,GLU,ASN,GLN".')
    parser.add_argument(
        '-mut_by_type',
        type=mutation_set_type,
        metavar='',
        help='Choose the residues to mutate by type.  Available options: '
        'polar, nonpolar, aromatic, neutral, basic, acidic, charged.')

    # Specify a custom destination type for each input residue:
    parser.add_argument(
        "-res_file",
        metavar='',
        dest="res_file",
        help='A file containing residue mutations. Each line defines the '
        'mutation targets for a single residue. The format of a '
        'line is "<chain>:<resnum><inscode(optional)> <comma-separated 3-letter '
        'mutations>", e.g. "A:24 ARG,GLU,ASN,GLN". This option overrides the "-mut" and '
        '"-residues" flags.')
    # Specify a custom list of mutations, each mutation can involve multile
    # residues:
    parser.add_argument(
        "-muts_file",
        metavar='',
        help='A file containing multiple residue mutations to perform. Each line '
        'defines ONE mutation, which may involve multiple residues and cannot combine '
        'with other mutations. Not compatible with the -concurrent and -sequential. '
        'The format of each line is a list of space-separated residue mutations: '
        '"<chain>:<resnum><inscode(optional)>" followed by "->" and then the '
        '3-letter residue name or 1-letter code, e.g. "A:24->SER A:25->R"'
        'It also supports loop insertion and deletions, e.g.:"'
        'A:24->+GLYA //insert 4 residues with the sequence GLYA after A:24'
        'A:24->-3 //delete 3 residues after A:24')

    # Define the structures of non-standard amino acids
    parser.add_argument(
        "-residue_structure",
        metavar='',
        default=None,
        help='A Mae file containing the non-standard residues that are mutation '
        'targets. The Mae file can have multiple structures. Each structure can be a '
        'free amino acid, or a amino acid capped by ACE and NME, and the structure must '
        'have the same 3-letter residue name as in the mutation definition.')
    # Alternatively, directly provide the rotamer library and associate Mae files
    # TODO

    parser.add_argument(
        '-calc',
        type=calculation_type,
        default=[],
        help='A list of properties to calculate for reference structure '
        'and mutant structure. Multiple calculations are separated by '
        'commas, e.g: "-calc pka,rotatable". NOTE: "prime_energy" will '
        'always be calculated when a Prime refinement is used, and '
        '"e_pot" will always be calculated when "python_minimize" refinement '
        'is used ')
    parser.add_argument(
        '-concurrent',
        type=int,
        metavar='',
        default=1,
        help='Maximum concurrent residue changes (default 1)')

    # NOTE This option is deprecated:
    parser.add_argument(
        "-sequential",
        action="store_true",
        default=False,
        help="Make concurrent changes to only sequential (neighbor) residues.")

    parser.add_argument(
        '-refine_mut',
        choices=refine_choices,
        default=protein.Refiner.PRIME_RESIDUE,
        help='Refinement method. Available options:'
        '"python_minimize" - Minimize using the Python Minimizer (ffld); '
        '"prime_residue" (default) - Refine with Prime side-chain prediction followed by residue minimization; '
        '"prime_minimize" (obsolete, same as prime_residue); '
        '"prime_sidechain" - Refine with Prime side-chain prediction; '
        '"prime_sidechain_cbeta" - Refine with Prime side-chain prediction with CA-CB vector sampling; '
        '"prime_sidechain_bb" - Refine with Prime side-chain prediction with backbone sampling.'
        '"prime_loop_prediction" - Refine with Prime loop prediction.')
    parser.add_argument(
        '-dist',
        '-d',
        type=float,
        default=0,
        metavar='Ang.',
        help='Cutoff distance from the mutated residues/ligand for which other '
        'residues should be included in refinement. The distance is measured '
        'from a reference arginine residue placed at the mutated position. '
        'Any residue containing an atom within the cutoff distance will be '
        'included in the refinement. Default: 0.0')
    parser.add_argument(
        '-loop',
        type=residue_type,
        default=None,
        help=
        'Two residues to define the loop. Each residue should be in the form '
        '<chain>:<resnum>. If there is no chain use "_". Example: "-loop A:6,A:10".'
    )
    parser.add_argument(
        '-loop2',
        type=residue_type,
        default=None,
        help=
        'Two residues to define the second loop in cooperative loop sampling.'
        'Each residue should be in the form '
        '<chain>:<resnum>. If there is no chain use "_". Example: "-loop2 A:22,A:27".'
    )
    parser.add_argument(
        '-loop_options',
        type=residue_type,
        default=None,
        help=
        'Options for loop refinement. Keyword/Value are separated by comma,e.g.:'
        'MAX_CA_MOVEMENT/4.0,RES_SPHERE/7.5,HOST/localhost:8,'
        'PROTOCOL/[LOOP_BLD][EXTENDED][LONG_LOOP_2]')

    # BIOLUM-150 (not implemented yet)
    '''
    parser.add_argument(
        '-always_refine', default='',
        help='ASL representing residues that should always be refined along with '
        'the mutated residues (do not include atom numbers in the ASL).'
    )
    parser.add_argument(
        '-never_refine', default='',
        help='ASL representing residues that should never be refined, even if they '
        'are within the cutoff distance from the mutated residues (do not include '
        'atom numbers in the ASL). NOTE: The mutated residues will always be '
        'refined, even if they match this ASL.'
    )
    '''

    parser.add_argument(
        '-ligand_asl',
        metavar='',
        help='ASL used to split the receptor and ligand for bound and unbound '
        'calculations. This is only needed if the "struct_file" is a single '
        'CT, receptor/ligand complex. If this argument is not supplied there '
        'will not be any bound/unbound calculations.  Note: The "-receptor_asl" '
        'cannot be used with this flag.')
    parser.add_argument(
        '-receptor_asl',
        metavar='',
        help='ASL used to split the receptor and ligand for bound and unbound '
        'calculations. This also triggers a binding affinity calculation, as '
        'with "-ligand_asl". Note: The "-ligand_asl" cannot be used with this '
        'flag.')
    parser.add_argument(
        "-not_idealized",
        action="store_true",
        dest="not_ideal",
        help='Turn off ideal reference mutations. Ideal reference mutations '
        'take the same methods which idealize a residue\'s bond length and '
        'angles in the mutant residue, and apply it to the reference. Turning '
        'this off will use the input bond length and angles in the reference '
        'when calculating the deltas.')
    # FIXME change the "not_ideal" to "idealize" (and action to "store_false")

    parser.add_argument(
        '-solvent',
        default='vsgb2.0',
        choices=['vsgb2.0', 'vsgb2.1', 'sgbnp', 'vac'],
        help='The sovlent model for Prime jobs. Available options are '
        'vsgb2.0 (default), vsgb2.1, sgbnp and vac.')

    parser.add_argument(
        "-use_membrane",
        action="store_true",
        dest="use_membrane",
        default=False,
        help='Use the implicit membrane, which has to be set up and saved '
        'in the input structure already.')

    # NOTE: This option is obsolete:
    parser.add_argument(
        "-fast",
        action="store_true",
        dest="fast",
        default=True,
        help='Use the fast Prime residue mutation backend (the default).')

    parser.add_argument(
        "-classic",
        action="store_true",
        dest="classic",
        default=False,
        help='Use the Python driver based residue mutation backend.')

    parser.add_argument(
        "-rigidmove",
        action="store_true",
        dest="rigidmove",
        default=False,
        help='Allow rigidbody movement between ligand and receptor.')

    parser.add_argument(
        "-refine_mutable_only",
        action="store_true",
        dest="refine_mutable_only",
        default=False,
        help='Refine only the mutable residues, i.e. all residues within'
        ' certain cutoff but not mutable will be ignored.')

    parser.add_argument(
        "-refine_all_mutable",
        action="store_true",
        dest="refine_all_mutable",
        default=False,
        help='Refine all mutable residues for each MC search step.')

    parser.add_argument(
        "-refine_unbound",
        action="store_true",
        dest="refine_unbound",
        default=False,
        help='Refine the unbound state using the specified refinement method.'
        'By default, the unbound state is defined as separating the ligand and '
        'the receptor from the refined bound state without further optimization.'
    )

    parser.add_argument(
        "-MC",
        action="store_true",
        dest="mc",
        default=False,
        help='Run Monte Carlo search with Prime residue mutation backend.')
    parser.add_argument(
        '-mc_nstep',
        dest="mc_nstep",
        type=int,
        default=2000,
        help='Number of Marte Carlo steps (default 2000).')
    parser.add_argument(
        '-mc_noutput',
        dest="mc_noutput",
        type=int,
        default=100,
        help='Number of top mutations being reported (default 100).')
    parser.add_argument(
        '-mc_temp',
        dest="mc_temp",
        type=float,
        default=298.0,
        help='The temperature for controlling acceptance rate (default 298).')
    parser.add_argument(
        "-mc_write_trj",
        action="store_true",
        dest="mc_write_trj",
        default=False,
        help=
        'Write MC sampling trajectory files. This option should be used with '
        '-no_cleanup to keep the files in the work directory.')
    parser.add_argument(
        '-mc_nterm',
        dest="mc_nterm",
        type=int,
        default=1000000,
        help='Terminate if no new move accepted in X moves (default X=1000000).'
    )
    parser.add_argument(
        '-mc_stability_weight',
        dest="mc_stability_weight",
        type=float,
        default=0.0,
        help='The weight of stability in the optimization score in the affinity '
        'calculation. The default is 0.0.')
    parser.add_argument(
        '-mc_affinity_weight',
        dest="mc_affinity_weight",
        type=float,
        default=1.0,
        help='The weight of affinity in the optimization score in the affinity '
        'calculation. The default is 1.0.')
    parser.add_argument(
        '-mc_stability_cutoff',
        dest="mc_stability_cutoff",
        type=float,
        default=30.0,
        help=
        'The cutoff to reject a mutation if the stability increases more than '
        ' a threshold. The default is 30.0.')
    parser.add_argument(
        '-mc_affinity_cutoff',
        dest="mc_affinity_cutoff",
        type=float,
        default=30.0,
        help=
        'The cutoff to reject a mutation if the affinity increases more than '
        ' a threshold. The default is 30.0.')

    parser.add_argument(
        '-mc_random_seed',
        dest="mc_random_seed",
        type=int,
        default=0,
        help='The random seed for MC job (default 0)')
    parser.add_argument(
        "-mc_random_start",
        action="store_true",
        dest="mc_random_start",
        default=False,
        help='Start the the Monte Carlo search from a randomly mutated structure'
    )

    parser.add_argument(
        '-minimizer_nbcutoff',
        type=float,
        metavar='',
        default=None,
        help='Experimental option. Sets the non-bonded cutoff for the Python '
        'Minimizer (default is 14A if Python Minimizer is used for refinement, '
        'and 40A if Prime is used for refinement).')

    parser.add_argument(
        "-no_cleanup",
        action="store_false",
        dest="cleanup",
        help="Do not clean up intermediate files.")
    parser.add_argument(
        "-NOJOBID",
        action="store_true",
        help="Do not run the script under job control")
    parser.add_argument(
        "-HOST",
        default="localhost",
        metavar="<host_list>",
        help="Run job remotely on the indicated host entry. (if -NJOBS is used, "
        "subjobs are run on the -subhost list)")
    parser.add_argument(
        "-subhost",
        default="",
        metavar="<subhost_list>",
        help="Run subjobs on this host. By default, same as the -HOST list.")
    parser.add_argument(
        "-NJOBS",
        default=1,
        type=int,
        metavar="<#jobs>",
        help=
        "Number of subjobs to generate. By default, the workflow is run as one job."
    )
    parser.add_argument(
        "-LOCAL",
        action="store_true",
        help=
        "Do not use a temporary directory. Keep files in the current directory."
    )
    parser.add_argument(
        "-WAIT",
        action="store_true",
        help="Do not return a prompt until the job completes.")
    # Handled by the top-level script:
    parser.add_argument(
        "-SAVE",
        action="store_true",
        help="Return zip archive of job directory at job completion.")

    # Put residue scanning in FEP preparation mode.  This will
    # return a Residue scanning score and structure designed to
    # be used by FEP ligand selectivity analysis.  It's help is hidden
    # because it should only be called as part of the appropriate FEP
    # workflow.
    parser.add_argument(
        "-fep_preparation", action="store_true", help=argparse.SUPPRESS)

    return parser


def residue_type(residues):
    """
    Return a list of residues (chain:residue_#).
    This method validates the passed in -residues using the Mutator method
    for this.

    :type residues: str
    :param residues: comma-separated string of chain:residue_#
    """

    return residues.split(',')


def mutation_type(residues):
    """
    Return a list of residues that will be mutated based on residue type.

    :type residues: str
    :param residues: comma-separated string of 3 letter residue names
    """
    residues = residues.split(',')

    # disable the checking to support nonstandard residues
    # try:
    #    protein.Mutator.validate_mutated_residues(residues)
    # except ValueError, msg:
    #    raise argparse.ArgumentTypeError(str(msg))

    return residues


def calculation_type(calculations):

    valid_calcs = protein.PropertyCalculator.AGGREGATE_CALCULATIONS
    calcs = calculations.split(',')
    for c in calcs:
        if not c in valid_calcs:
            msg = 'Invalid calculation: %s\n' % c
            msg += 'Allowed types: %s' % ', '.join(valid_calcs)
            raise argparse.ArgumentTypeError(msg)

    return calcs


def mutation_set_type(mut_type):
    """
    Return a list of residues that will be mutated based on residue type.

    """
    valid_mut_types = [
        'all', 'polar', 'nonpolar', 'aromatic', 'neutral', 'charged', 'basic',
        'acidic'
    ]

    if not mut_type in valid_mut_types:
        msg = 'Invalid mutation by type: %s\n' % mut_type
        msg += 'Allowed types: %s' % ', '.join(valid_mut_types)
        raise argparse.ArgumentTypeError(msg)

    if mut_type == 'all':
        mutations = protein.ALL_RESIDUES
    elif mut_type == 'polar':
        mutations = protein.POLAR_RESIDUES
    elif mut_type == 'nonpolar':
        mutations = protein.NONPOLAR_RESIDUES
    elif mut_type == 'aromatic':
        mutations = protein.AROMATIC_RESIDUES
    elif mut_type == 'neutral':
        mutations = protein.NEUTRAL_RESIDUES
    elif mut_type == 'charged':
        mutations = protein.CHARGED_RESIDUES
    elif mut_type == 'basic':
        mutations = protein.BASIC_RESIDUES
    elif mut_type == 'acidic':
        mutations = protein.ACIDIC_RESIDUES

    return mutations


def get_prime_mutations_list(file_name):
    """
    This function reads the output maegz file from Prime, extracts
    and returns a list of the mutations (stored in the property
    's_bioluminate_Mutations') in the order they occur in the Prime
    output maegz file

    :param file_name: name of file from which mutations to be extracted
    :type file_name: string
    :returns: list of tuples indicating mutation information
    :rtype: list of tuples

    """
    regex = PRIME_MUTANT_RE
    prime_mutations_list = []
    for st in structure.StructureReader(file_name):
        try:
            desc = st.property['s_bioluminate_Mutations']
        except KeyError:
            continue
        changes = []
        if desc == 'NONE' or desc == 'MINIMIZED':
            continue
        for change_str in desc.split(','):
            match = regex.search(change_str)
            if match:
                chain = match.group('chain')
                resnum = match.group('resnum')
                inscode = match.group('inscode') or ' '
                new_resname = match.group('new_resname')
                changes.append((chain, int(resnum), inscode, new_resname))
            else:
                raise RuntimeError(
                    "Invalid mutation recognized: %s" % change_str)
        prime_mutations_list.append(changes)

    if not prime_mutations_list:
        raise RuntimeError("No mutations found in the specified file")

    return prime_mutations_list


def write_reordered_mutants(file_name, muts_file):
    """
    This method compares the list of mutations requested by the user (as
    in muts_file and reorders the cts in the output maegz file and
    replaces it with a new output maegz file

    :param file_name: name of maegz file with mutants
    :type file_name: string

    """
    prime_mutations_list = get_prime_mutations_list(file_name)
    mutations_list = protein.Mutator.convert_muts_file(
        muts_file, regex=SINGLE_MUTATION_RE)

    if len(prime_mutations_list) != len(mutations_list):
        info('== NOTE: Difference in # of requested & generated mutants ==')

    unordered_outfile = file_name
    reordered_outfile = 'ordered_output.maegz'

    #Read the unordered mutant file into a list
    unordered_sts = []
    num_unmut = 0
    for st in structure.StructureReader(unordered_outfile):
        unordered_sts.append(st)
        if st.property['s_bioluminate_Mutations'] == 'NONE' or\
             st.property['s_bioluminate_Mutations'] == 'MINIMIZED':
            num_unmut += 1

    # Write the mutants file in correct order
    with structure.StructureWriter(reordered_outfile) as writer:
        for num in range(num_unmut):
            writer.append(unordered_sts[num])
        for mut_list in mutations_list:
            try:
                prime_mut_index = prime_mutations_list.index(mut_list)
            except ValueError:
                info('Mutation %s has not been generated' % mut_list)
            else:
                writer.append(unordered_sts[num_unmut + prime_mut_index])

    shutil.copy(reordered_outfile, unordered_outfile)
    os.remove(reordered_outfile)


def get_parser_args(parser):
    """
    Get the args from a parser intended for use with residue_scanning_backend.
    This will check to make sure the struct file is available.

    """
    try:
        args = parser.parse_args()
    except IOError as msg:
        parser.error(str(msg))

    debug("CWD: %s" % os.getcwd())

    if not os.path.isfile(args.structures):
        msg = 'Structure file not found: "%s"' % args.structures
        parser.error(msg)

    if args.res_file and not os.path.isfile(args.res_file):
        msg = 'Residue parameter file not found: "%s"' % args.res_file
        parser.error(msg)

    if args.muts_file and not os.path.isfile(args.muts_file):
        msg = 'Mutations parameter file not found: "%s"' % args.muts_file
        parser.error(msg)

    if not args.residues and not args.res_file and not args.muts_file:
        msg = 'One of the following must be specified: -residues, -res_file, -muts_file.'
        parser.error(msg)

    if args.residues:
        if not args.mut and not args.mut_by_type:
            msg = 'One of the following must be specified: -mut, -mut_by_type.'
            parser.error(msg)


#   this limit should be lifted because MC search can do arbitrary number of mutations
#    if args.concurrent > 99:
#        msg = "-concurrent: Up to 99 maximum concurrent mutations are supported."
#        parser.error(msg)

    if args.mc and args.concurrent == 1:
        msg = "concurrent must be set to greater than 1 for MC job."
        parser.error(msg)

    if not args.mc and (args.concurrent > 1 and args.muts_file):
        msg = "-concurrent option is not supported with -muts_file for non-MC jobs."
        parser.error(msg)

    if args.mc_write_trj and args.cleanup:
        msg = "-mc_write_trj should be used with -no_cleanup."
        parser.error(msg)

    if args.ligand_asl and args.receptor_asl:
        msg = 'The "-ligand_asl" and "-receptor_asl" cannot be used together.'
        parser.error(msg)

    if args.refine_mut != protein.Refiner.PYTHON_MINIMIZE and \
                not (os.path.isfile(PRIME) or os.path.isfile(PRIMEEXE)):
        msg = "Could not find $SCHRODINGER/prime; make sure Prime is installed"
        parser.error(msg)

    if args.refine_mut == protein.Refiner.PYTHON_MINIMIZE:
        args.fast = False
        args.classic = True

    if args.loop and args.refine_mut != protein.Refiner.PRIME_LOOP_PRED:
        msg = 'Refinement method must be Prime loop prediction when -loop is specified'
        parser.error(msg)

    if args.loop:
        args.classic = True

    if args.fast and args.classic:
        args.fast = False

    if args.mc and args.classic:
        msg = '-MC and -classic are not compatible.'
        parser.error(msg)

    if args.residue_structure and not (args.fast or args.mc):
        msg = 'Non-standard residues are only supported with fast mode.'
        parser.error(msg)

    return args


def submit_under_jobcontrol(args, script_args):
    """
    Submit this script under job control.

    """
    # Launch the backend under job control:

    prog_name = "Residue Mutation"

    scriptlauncher = launcher.Launcher(
        script=BACKEND,
        jobname=args.jobname,
        runtoplevel=True,  # FIXME should be False
        # FIXME It just does not seem right to have runtoplevel
        # be set to True. Does this really fix remote jobs?????
        # It may have negitive side-effects. The GUI runs the
        # backend via $SCHRODINGER/run anyway.
        prog=prog_name,
        wait=args.WAIT,
        no_redirect=False,
    )
    # NOTE: -SAVE is passed from toplevel to jlaunch via
    # SCHRODINGER_SAVE_JOBDIR.

    # Add input file to jlaunch input file list:
    scriptlauncher.addInputFile(args.structures)

    # NOTE: Output file will be added at run-time

    if args.res_file:
        scriptlauncher.addInputFile(args.res_file)
    if args.muts_file:
        scriptlauncher.addInputFile(args.muts_file)
    if args.residue_structure:
        scriptlauncher.addInputFile(args.residue_structure)

    # Add script arguments:
    # Since the backend will already be running under job control.
    script_args.append("-NOJOBID")

    # NOTE: We should NOT remove the -LOCAL argument from the
    # args list, because the backend will need to use it.

    # Tell the driver where to run the subjobs (if the user did not specify
    # -subhost):
    if args.NJOBS > 1 and not args.subhost:
        host_list = jobcontrol.get_command_line_host_list()
        host_string = jobcontrol.host_list_to_str(host_list)
        script_args += ["-subhost", host_string]

    scriptlauncher.addScriptArgs(script_args)

    scriptlauncher.launch()


def _run_backend(args, backend=None):
    """
    This function is run after a BIOLUMINATE_MAIN licence is checked out
    It is designed to either run a single processor job
    (args.subjob = False) or run a subjob of a multiprocessor job
    (args.subjob = True).  In either case license checking should be done
    before this subroutine is called.

    """
    start_time = time.time()

    if backend:
        backend.setJobProgress(description="Generating mutations")

    info(
        'Backend command (used for starting the backend in the job directory):')
    info('%s\n' % subprocess.list2cmdline(sys.argv))

    struct_file = jobcontrol.get_runtime_path(args.structures)

    # FIXME: Right now we are assuming that the input is a single
    # CT, protein/ligand complex. Allow for multiple CTs.
    original_st = structure.Structure.read(struct_file)

    info("Job name: %s" % args.jobname)

    if args.ligand_asl or args.receptor_asl:
        info('This job will run binding affinity calculations')
    else:
        info('This job will run protein stability calculations')

    info("Refinement method: %s" % args.refine_mut)
    info("Will use refinement distance of %fA" % args.dist)

    min_nbcutoff = args.minimizer_nbcutoff
    if min_nbcutoff is None:
        if args.refine_mut == protein.Refiner.PYTHON_MINIMIZE:
            min_nbcutoff = 14.0
        else:
            # See BIOLUM-1435
            min_nbcutoff = 40.0
    info(
        "Will use non-bonded cutoff for Python Minimizer of %fA" % min_nbcutoff)

    res_file = args.res_file
    if res_file:
        res_file = jobcontrol.get_runtime_path(res_file)
        info("Runtime residue mutations file path: %s" % res_file)

    muts_file = args.muts_file
    if muts_file:
        muts_file = jobcontrol.get_runtime_path(muts_file)
        info("Runtime mutations combinations file path: %s" % muts_file)

    info("Analyzing requested mutations list...")
    info("  Maximum number of concurrent mutations: %i" % args.concurrent)

    idealize = (not args.not_ideal)

    if args.muts_file:
        try:
            mutations_list = protein.Mutator.convert_muts_file(
                muts_file, regex=SINGLE_MUTATION_RE)
        except RuntimeError as msg:
            error(str(msg))
            sys.exit(1)
        mutator = protein.Mutator(original_st, [], args.concurrent,
                                  args.sequential, idealize)
    else:
        try:
            if res_file:
                requested_mutations = protein.Mutator.convert_res_file(
                    res_file, regex=MULTI_MUTATION_RE)
            else:
                if args.mut:
                    muts = args.mut
                else:
                    muts = args.mut_by_type
                requested_mutations = protein.Mutator.convert_residue_list(
                    args.residues, muts, regex=MULTI_MUTATION_RE)
        except RuntimeError as msg:
            parser.error(str(msg))

        concurrent = args.concurrent
        # MC job will handle concurrent itself
        if args.mc:
            concurrent = 1
        mutator = protein.Mutator(original_st, requested_mutations, concurrent,
                                  args.sequential, idealize)
        mutations_list = mutator.calculateMutationsList()

    # filter out crosslinked residues, BIOLUM-2632
    templist = []
    skipped_res = []
    for mut_changes in mutations_list:
        skip = False
        for chain, resnum, inscode, new_resname in mut_changes:
            resid = "%s:%s%s" % (chain, resnum, inscode)
            if resid in skipped_res:
                skip = True
                break
            if psp_util.is_residue_crosslinked(original_st, resid) == 2:
                info("  Warning: crosslinked residue %s will be skipped." %
                     resid)
                skip = True
                skipped_res.append(resid)
                break
        if not skip:
            templist.append(mut_changes)
    mutations_list = templist

    if len(mutations_list) == 0:
        error("The mutation combinations list is empty")
        sys.exit(1)

    # Print all mutation combinations:
    info("  List of %i potential mutants" % len(mutations_list))
    for i, mut_changes in enumerate(mutations_list, start=1):
        mut_strings = []
        for chain, resnum, inscode, new_resname in mut_changes:
            if inscode == " ":
                inscode = ""
            mut_str = "%s:%s%s->%s" % (chain, resnum, inscode, new_resname)
            mut_strings.append(mut_str)
        info("  %5i %s" % (i, ", ".join(mut_strings)))

    calculations = args.calc
    # BIOLUM-1276 Even if the user did not specify it, we need to calculate
    # the minimizer energy in order to calculate the affinity:
    # BIOLUM-2981: only do it when no property is specified in -classic mode
    if args.classic:
        if not 'e_pot' in calculations:
            calculations.append("e_pot")
    # We need to pull out the prime_energy calculation since it takes so long.
    # This is later run using the self.calculate_prime_energy property created
    # in the __init__.
    if 'prime_energy' in calculations:
        calculations.remove('prime_energy')

    out_file = '%s-out.maegz' % args.jobname  # refined mutant
    if os.path.exists(out_file):
        os.remove(out_file)
    # refined reference (self-mutated) BIOLUM-1216
    refs_file = '%s-refs.maegz' % args.jobname
    if os.path.exists(refs_file):
        os.remove(refs_file)

    # Decide whether to do Prime residue mutation or not
    # If yes, do it and return
    if args.fast or args.mc \
           and (args.refine_mut != 'prime_loop_prediction'):
        # does not support loop prediction yet
        info("\nRun Prime Residue Mutation...")
        ligand_asl = None
        if args.ligand_asl:
            ligand_asl = args.ligand_asl
        elif args.receptor_asl:
            ligand_asl = 'NOT (%s)' % args.receptor_asl

        # Don't calculate the surface complementarity when there is no ligand
        # (or receptor) defined:
        if 'vdw_surf_comp' in calculations and not ligand_asl:
            calculations.remove('vdw_surf_comp')

        prs = PrimeResidueScan(
            original_st,
            mutations_list,
            args.jobname,
            nojobid=True,  # Run scans without (extra) job control
            refine_mut=args.refine_mut,
            ligand_asl=ligand_asl,
            residue_structure=args.residue_structure,
            dist=args.dist,
            mutable_only=args.refine_mutable_only,
            all_mutable=args.refine_all_mutable,
            rigidmove=args.rigidmove,
            refine_unbound=args.refine_unbound,
            solvent=args.solvent,
            use_membrane=args.use_membrane,
            mc_job=args.mc,
            nstep=args.mc_nstep,
            noutput=args.mc_noutput,
            temperature=args.mc_temp,
            write_trj=args.mc_write_trj,
            nterm=args.mc_nterm,
            stab_weight=args.mc_stability_weight,
            aff_weight=args.mc_affinity_weight,
            stab_cutoff=args.mc_stability_cutoff,
            aff_cutoff=args.mc_affinity_cutoff,
            concurrent=args.concurrent,
            random_seed=args.mc_random_seed,
            random_start=args.mc_random_start,
            cleanup=args.cleanup,
            fep_preparation=args.fep_preparation)
        t1 = time.time()
        prs.run()
        t2 = time.time()
        info('TIME for Prime Residue Mutation %f seconds' % (t2 - t1))

        # take Prime output and calculate other properties
        if len(prs.outfiles) != 2 or \
               (not os.path.exists(prs.outfiles[0])) or \
               (not os.path.exists(prs.outfiles[1])):
            error("Missing outputs from Prime Residue Mutation!")
            sys.exit(1)

        prs.calcDeltaProp(calculations, min_nbcutoff, refs_file, out_file)

        if args.cleanup:
            prs.clean()

        if backend:
            # If running under job control
            backend.addOutputFile(out_file)
            backend.addOutputFile(refs_file)
            backend.setStructureOutputFile(out_file)

        end_time = time.time()
        info('\nElapsed time: %f seconds.' % (end_time - start_time))
        return

    # Now continue to discrete mode

    # Create the writer so we can append all new structures
    writer = structure.StructureWriter(out_file)
    refs_writer = structure.StructureWriter(refs_file)

    writer.append(original_st)

    # The arg ,ligand_asl/receptor_asl, will let us know there is a
    # bound/unbound calc we need to do.
    if args.ligand_asl:
        scan_type = Scanner.LIGAND_SUBUNIT
        subunit_asl = args.ligand_asl
    elif args.receptor_asl:
        scan_type = Scanner.PROTEIN_SUBUNIT
        subunit_asl = args.receptor_asl
    if args.ligand_asl or args.receptor_asl:
        scanner = Scanner(
            original_st,
            args.refine_mut,
            scan_type,
            subunit_asl,
            solvent=args.solvent,
            use_membrane=args.use_membrane,
            loop=args.loop,
            loop_opt=args.loop_options,
            loop2=args.loop2,
            cleanup=args.cleanup,
            nbcutoff=min_nbcutoff,
        )
        for subunit in scanner.ref_subunits:
            if subunit:
                writer.append(subunit)
    else:
        scanner = Scanner(
            original_st,
            args.refine_mut,
            solvent=args.solvent,
            use_membrane=args.use_membrane,
            loop=args.loop,
            loop_opt=args.loop_options,
            loop2=args.loop2,
            cleanup=args.cleanup,
            nbcutoff=min_nbcutoff,
        )

    # Determine which residues would NOT be idealized:
    # (they are reported to the user at the end of the workflow)
    skipped_residues = []
    if idealize:
        input_residue_list = []
        for changes in mutations_list:
            for chain, resnum, inscode, new_resname in changes:
                res_id = (chain, resnum, inscode)
                if not res_id in input_residue_list:
                    input_residue_list.append(res_id)

        for chain, resnum, inscode in input_residue_list:
            # Get the residue that is going to be mutated
            residue = original_st.findResidue('%s:%s%s' % (chain, resnum,
                                                           inscode))
            resname = residue.pdbres.strip()

            # Catch any unsupported residues and do not try to mutate
            # them. Add this information to the log in the summary.
            if not resname in protein.Mutator.SUPPORTED_BUILD_RESIDUES:
                skipped_residues.append((chain, resnum, resname))
                continue

            # NOTE: The actual self-mutation will be performed in the Mutator class
            # See BIOLUM-1215

            # Set the total steps to the total mutations being run + 1 more
            # for the initial reference property calculations
    total_steps = len(mutations_list) + 1
    step = 1
    if backend:
        backend.setJobProgress(step, total_steps)

    info("\nGenerating mutations...")

    skipped_mutations = 0
    # Generate all of the mutations for the reference and if there
    # is a refinement add them to the refinement queue
    # Each mutation yielded by this function contains a mutated CT
    # which is raw; that is, it is unrefined in any way.

    for changes in mutations_list:
        # Each iteration of this loop represents a MUTATION. Each mutation
        # can consist or one or more residue changes.
        jobname = args.jobname + "_%i" % step

        change_strings = []
        for chain, resnum, inscode, new_resname in changes:
            if inscode == " ":
                inscode = ""
            change_strings.append("%s:%s%s->%s" % (chain, resnum, inscode,
                                                   new_resname))
        info("\nProcessing mutant: %s   %s" % (jobname,
                                               ", ".join(change_strings)))

        # This is where we refine the mutation (and the reference structure the same way),
        # and generate the properties:
        try:
            # BIOLUM-150 Add always_refine and never_refine here also:
            mutation = mutator.getMutationFromChanges(changes)
            ref_st, sts = scanner.calculateOneMutation(mutation, args.dist,
                                                       calculations, jobname)
            info("Writing structure(s) to file...")
            if ref_st:
                refs_writer.append(ref_st)
            for st in sts:
                writer.append(st)

        except Exception as msg:
            error("ERROR: %s" % msg)
            warning("Skipping mutation...")
            skipped_mutations += 1
        finally:
            # Increment the step
            step += 1
            info(str(100 * step / total_steps) + "% complete\n")
            if backend:
                backend.setJobProgress(step, total_steps)

    # Ev:120121 Make sure the progress gets to the end:
    if backend:
        backend.setJobProgress(step, total_steps)

    num_out = len(mutations_list) - skipped_mutations

    if num_out == 0:
        error("ERROR: No mutants were generated. Exiting.")
        sys.exit(1)

    info('\nNumber of mutants generated: %i' % num_out)

    if skipped_residues:
        warning('\nThe following reference residues could not be idealized.')
        warning('The residues were mutated but the comparisons were done ')
        warning('with a non-idealized reference.')
        for residue in skipped_residues:
            warning('    %s:%s %s' % residue)

    writer.close()
    refs_writer.close()

    info('\nOutput file: %s' % out_file)
    debug('Reference (self-mutated refined structures) file: %s' % refs_file)

    if backend:
        # If running under job control
        backend.addOutputFile(out_file)
        backend.addOutputFile(refs_file)
        backend.setStructureOutputFile(out_file)

    end_time = time.time()
    info('\nElapsed time: %f seconds.' % (end_time - start_time))


# Tampering with licensing is a violation of the license agreement


def check_license(args):
    """
    Check to see if the correct licenses are present for a given set
    if input settings.  If a license will need to be checked in
    then return the handle to do so, otherwise return None.
    This should be called at the start of every job and subjob.  If this
    is a subjob, then no license is checked out.  If this is a either
    a single processor job OR the driver for a multiprocessor job then
    a FEP_GPGPU token is checked for if the -fep_preparation flag
    was used or a BIOLUMINATE_MAIN token is checked out for the duration
    of the multi-processor job otherwise.
    """
    license = None
    # BIOLUM-1630 Subjob do not need to check-out a new license
    if args.subjob:
        pass
    elif args.fep_preparation:
        if not biolicense.is_license_available(biolicense.FEP_GPGPU):
            error("ERROR: An FEP license is required to run " +
                  "Residue Scanning in FEP Preparation Mode")
            sys.exit(1)
    else:
        license = biolicense.License(biolicense.BIOLUMINATE_MAIN)
        if not license.isValid():
            error("ERROR: No BIOLUMINATE_MAIN license token is available. "
                  "You will not be able to run the calculation.")
            sys.exit(1)
        # The license will be checked-in when the workflow finishes.
    return license


def run_backend(args):
    """
    Run the backend job.  The called subroutine _run_backend contains the
    actual job processing machinery, but the exception handling and
    licensing are handled here.
    (We are either already running under job control or -NOJOBID was used)
    """

    # Tampering with licensing is a violation of the license agreement
    license = check_license(args)
    backend = jobcontrol.get_backend()

    # Run the backend, and make sure the license is checked in, even if it
    # fails:
    try:
        _run_backend(args, backend)
    except:
        # Make sure we reset the progress bar if the job fails (Ev:128964)
        if backend:
            backend.setJobProgress()
        # FIXME it would be better to implement a universal implementation
        # similar to this in AppFramework itself.

        # We still need to show the traceback:
        raise
    finally:
        if license is not None:
            license.checkin()


def merge_subjob_results(main_jobname, sub_job_names, opt_affinity, fast_mode,
                         args):
    """
    Merge the results from the given subjobs. Returns the combined
    output file paths.
    """

    combined_outfile = main_jobname + "-out.maegz"
    combined_refsfile = main_jobname + "-refs.maegz"
    combined_writer = structure.StructureWriter(combined_outfile)
    combined_refs_writer = structure.StructureWriter(combined_refsfile)

    num_written_muts = 0
    mutation_sets = set()
    for subjob_num, subjobname in enumerate(sub_job_names, start=1):
        subjob_out_file = '%s-out.maegz' % subjobname
        subjob_refs_file = '%s-refs.maegz' % subjobname

        if not os.path.isfile(subjob_out_file):
            warning("WARNING: Subjob output file does not exist: %s" %
                    subjob_out_file)
            # It's OK, because we are now simply skipping failed subjobs
            continue

        if not os.path.isfile(subjob_refs_file):
            error("ERROR: Subjob references file does not exist: %s" %
                  subjob_refs_file)
            sys.exit(1)

        for i, st in enumerate(
                structure.StructureReader(subjob_out_file), start=1):
            if fast_mode:
                # fast mode always returns one structure per output.
                mutval = st.property['s_bioluminate_Mutations']
                # original starting structure
                if mutval == 'NONE' or mutval == 'MINIMIZED':
                    if subjob_num == 1:
                        combined_writer.append(st)
                    continue
                if mutval in mutation_sets:
                    continue
                mutation_sets.add(mutval)
                num_written_muts += 1
                index = st.title.rindex('_')
                st.title = st.title[0:index + 1] + str(num_written_muts)
                combined_writer.append(st)
            else:
                # in each subjob, the first output is the original starting structure.
                # stability calculation returns one structure per output, but affinity
                # calculation returns 3 structures per output.
                if opt_affinity:
                    is_ref_st = (i <= 3)
                else:
                    is_ref_st = (i == 1)
                if is_ref_st:
                    if subjob_num == 1:
                        combined_writer.append(st)
                else:
                    combined_writer.append(st)
                    # Report each mutation only once:
                    if st.property[TYPE_PROP] == COMPLEX_TYPE:
                        num_written_muts += 1

        for st in structure.StructureReader(subjob_refs_file):
            combined_refs_writer.append(st)

    combined_writer.close()
    combined_refs_writer.close()

    # sort the outputs from mc subjobs
    if args.mc and not args.fep_preparation:
        temp_outfile = main_jobname + "-out-temp.maegz"
        temp_writer = structure.StructureWriter(temp_outfile)
        props = []
        for i, st in enumerate(
                structure.MaestroTextReader(combined_outfile), start=1):
            if st.property['s_bioluminate_Mutations'] == 'NONE':
                prop = -2.0E10
            elif st.property['s_bioluminate_Mutations'] == 'MINIMIZED':
                prop = -1.0E10
            else:
                if opt_affinity:
                    prop = st.property[DELTA_AFFINITY_PROP]
                else:
                    prop = st.property[DELTA_STABILITY_PROP]
            props.append((i, prop))
        sorted_props = sorted(props, key=lambda x: x[1])
        for (i, prop) in sorted_props:
            st = structure.Structure.read(combined_outfile, index=i)
            temp_writer.append(st)
        temp_writer.close()
        shutil.copy(temp_outfile, combined_outfile)
        os.remove(temp_outfile)

    info("Number of combined mutants written out: %i" % num_written_muts)
    info("\nOutput file: %s" % combined_outfile)
    debug("Reference (self-mutated refined structures) file: %s" %
          combined_refsfile)

    return (combined_outfile, combined_refsfile)


def multijob_driver(args):
    """
    For a given set of inputs paramters, args, stored as {argparse.Namespace}
    prepare a multiprocessor job
    The actual task of launching subjobs is handled by _multijob_driver, but
    this will determine whether the jobs is running under jobcontrol
    and check out the applicable driver licenses before creating and running
    the subjobs.
    """
    backend = jobcontrol.get_backend()
    # Tampering with licensing is a violation of the license agreement
    license = check_license(args)
    _multijob_driver(args, backend)
    if license is not None:
        license.checkin()


def _multijob_driver(args, backend=None):
    """
    Split the workflow into the specified number of subjobs, submit each
    subjob under job control, wait until all subjobs complete, and combine
    the results.
    """

    script_args = sys.argv[1:]
    info('Main command: %s\n' % subprocess.list2cmdline(sys.argv))

    # Read the -HOST value, and convert to list (to pass to JobDJ):
    subhost_list = jobcontrol.host_str_to_list(args.subhost)
    info("Sub host list: %s" % subhost_list)
    jdj = queue.JobDJ(
        verbosity="normal", hosts=subhost_list, max_failures=queue.NOLIMIT)
    # FIXME currently restartability is not supported
    restarting = False

    if args.res_file:
        res_file = jobcontrol.get_runtime_path(args.res_file)
        info("Runtime residue mutation file path: %s" % res_file)

    if args.muts_file:
        muts_file = jobcontrol.get_runtime_path(args.muts_file)
        info("Runtime mutation combinations file path: %s" % muts_file)

    # Read the input structure:
    struct_file = jobcontrol.get_runtime_path(args.structures)
    original_st = structure.Structure.read(struct_file)

    # Make sure the mutations are specified correctly:
    info("Analyzing requested mutations list...")
    info("  Maximum number of concurrent mutations: %i" % args.concurrent)
    if args.muts_file:
        try:
            mutations_list = protein.Mutator.convert_muts_file(
                muts_file, regex=SINGLE_MUTATION_RE)
        except RuntimeError as msg:
            error(str(msg))
            sys.exit(1)
    else:
        try:
            if args.res_file:
                requested_mutations = protein.Mutator.convert_res_file(
                    res_file, regex=MULTI_MUTATION_RE)
            else:
                if args.mut:
                    muts = args.mut
                else:
                    muts = args.mut_by_type
                requested_mutations = protein.Mutator.convert_residue_list(
                    args.residues, muts, regex=MULTI_MUTATION_RE)
        except RuntimeError as msg:
            error(str(msg))
            sys.exit(1)

        # Generate all mutation combinations (BIOLUM-1275):
        idealize = (not args.not_ideal)
        concurrent = args.concurrent
        # MC job will handle concurrent itself
        if args.mc:
            concurrent = 1
            info(
                "  MC search will handle concurrent mutations later in Prime backend"
            )

        mutator = protein.Mutator(original_st, requested_mutations, concurrent,
                                  args.sequential, idealize)
        mutations_list = mutator.calculateMutationsList()

    # filter out crosslinked residues, BIOLUM-2632, BIOLUM-2749
    templist = []
    skipped_res = []
    for mut_changes in mutations_list:
        skip = False
        for chain, resnum, inscode, new_resname in mut_changes:
            resid = "%s:%s%s" % (chain, resnum, inscode)
            if resid in skipped_res:
                skip = True
                break
            if psp_util.is_residue_crosslinked(original_st, resid) == 2:
                info("  Warning: crosslinked residue %s will be skipped." %
                     resid)
                skip = True
                skipped_res.append(resid)
                break
        if not skip:
            templist.append(mut_changes)
    mutations_list = templist

    if len(mutations_list) == 0:
        error("ERROR: The mutant list is empty")
        sys.exit(1)

    # Print all mutation combinations:
    info("  List of %i potential mutants" % len(mutations_list))
    for i, mut_changes in enumerate(mutations_list, start=1):
        mut_strings = []
        for chain, resnum, inscode, new_resname in mut_changes:
            if inscode == " ":
                inscode = ""
            mut_str = "%s:%s%s->%s" % (chain, resnum, inscode, new_resname)
            mut_strings.append(mut_str)
        info("  %5i %s" % (i, ", ".join(mut_strings)))

    info("")

    # Determine how many subjobs to split into:
    numjobs = min(len(mutations_list), args.NJOBS)
    if numjobs < args.NJOBS:
        warning(
            "Creating fewer number of jobs: %i (the number of residues to mutate)"
            % numjobs)
    else:
        info("Creating %i subjobs" % numjobs)

    # For MC, each subjob will start with differnt random seed, do a fraction
    # of requested MC steps, and return a fraction of total outputs.
    if args.mc:
        random_seeds = [x for x in range(numjobs)]
        num_step_per_job = calc_items_per_job(args.mc_nstep, numjobs)
        nsteps = list(block_size(args.mc_nstep, num_step_per_job))
        num_output_per_job = calc_items_per_job(args.mc_noutput, numjobs)
        noutputs = list(block_size(args.mc_noutput, num_output_per_job))
    else:
        # some dummy numbers for non MC job
        random_seeds = [x for x in range(numjobs)]
        nsteps = [x for x in range(numjobs)]
        noutputs = [x for x in range(numjobs)]

    if args.mc:
        mutations_by_job = [mutations_list for x in range(numjobs)]
    else:
        # Determine which subjob will mutate which residues:
        num_muts_per_job = calc_items_per_job(len(mutations_list), numjobs)
        mutations_by_job = list(chunks(mutations_list, num_muts_per_job))

    # NOTE: This algorithm may produce less subjobs than requested.
    # For example, 5 mutations will be split up into 3 subjobs (2, 2, 1).

    main_jobname = args.jobname
    restart_file = main_jobname + '.restart'

    # Read the structures from the input file.
    sub_job_names = []

    for i, (seed, nstep, noutput, job_mutations) in enumerate(
            zip(random_seeds, nsteps, noutputs, mutations_by_job), start=1):
        sub_job_name = main_jobname + "-%03d" % (i)  # fill with zeros
        sub_job_names.append(sub_job_name)
        if restarting:
            continue

        info("Preparing subjob %i..." % i)

        run = os.path.join(os.environ["SCHRODINGER"], "run")
        backend_exe = "residue_scanning_backend.py"

        cmd = [run, backend_exe]

        cmd += ["-jobname", sub_job_name, '-subjob', str(i)]

        # Use the same input (pre-mutated) structure
        cmd.append(args.structures)

        if not args.mc:
            info("Number of mutants that will be created by this subjob: %i" %
                 len(job_mutations))

        # Write a mutations file for this subjob:
        mut_file = sub_job_name + "-muts.txt"
        fh = open(mut_file, "w")
        info("Mutants for this subjob:")
        for mut_changes in job_mutations:
            # Write all residues involved in this mutation as a row to the
            # subjob input file:
            change_strings = []
            for chain, resnum, inscode, new_resname in mut_changes:
                if inscode == " ":
                    inscode = ""
                change_strings.append("%s:%s%s->%s" % (chain, resnum, inscode,
                                                       new_resname))
            line = " ".join(change_strings)
            fh.write(line + "\n")
            info("  %s" % ", ".join(change_strings))
        fh.close()
        info("")

        cmd += ["-muts_file", mut_file]

        if args.residue_structure:
            cmd += ["-residue_structure", args.residue_structure]

        if args.calc:
            cmd += ["-calc", ",".join(args.calc)]

        # NOTE: The -concurrent and -sequential flags are not used here on purpose
        # but MC job needs -concurrent.
        if args.mc:
            cmd.append("-MC")
            cmd += ["-concurrent", str(args.concurrent)]
            cmd += ["-mc_nstep", str(nstep)]
            cmd += ["-mc_noutput", str(noutput)]
            cmd += ["-mc_random_seed", str(seed)]
            cmd += ["-mc_temp", str(args.mc_temp)]
            cmd += ["-mc_stability_weight", str(args.mc_stability_weight)]
            cmd += ["-mc_affinity_weight", str(args.mc_affinity_weight)]
            cmd += ["-mc_stability_cutoff", str(args.mc_stability_cutoff)]
            cmd += ["-mc_affinity_cutoff", str(args.mc_affinity_cutoff)]
            cmd += ["-mc_nterm", str(args.mc_nterm)]
            if args.mc_write_trj:
                cmd.append("-mc_write_trj")
            if args.mc_random_start:
                cmd.append("-mc_random_start")

        # both MC and fast job use the same output
        fast_mode = True
        if args.classic:
            cmd.append("-classic")
            fast_mode = False

        cmd += ["-refine_mut", args.refine_mut]
        cmd += ["-dist", str(args.dist)]
        cmd += ["-solvent", str(args.solvent)]
        if args.loop:
            cmd += ["-loop", ','.join(args.loop)]
            if args.loop_options:
                cmd += ["-loop_options", ','.join(args.loop_options)]
            if args.loop2:
                cmd += ["-loop2", ','.join(args.loop2)]
        if args.refine_mutable_only:
            cmd += ["-refine_mutable_only"]
        if args.refine_all_mutable:
            cmd += ["-refine_all_mutable"]
        if args.rigidmove:
            cmd += ["-rigidmove"]
        if args.use_membrane:
            cmd += ["-use_membrane"]
        if args.refine_unbound:
            cmd += ["-refine_unbound"]

        # BIOLUM-150
        # if args.always_refine:
        #    cmd += ["-always_refine", args.always_refine]
        # if args.never_refine:
        #    cmd += ["-never_refine", args.never_refine]

        if args.minimizer_nbcutoff is not None:
            cmd += ["-minimizer_nbcutoff", str(args.minimizer_nbcutoff)]

        if args.ligand_asl:
            cmd += ["-ligand_asl", args.ligand_asl]
        if args.receptor_asl:
            cmd += ["-receptor_asl", args.receptor_asl]

        # Used to determine how to merge subjob results:
        opt_affinity = args.ligand_asl or args.receptor_asl
        if args.mc_stability_weight >= 0.99 and args.mc_affinity_weight <= 0.01:
            opt_affinity = False

        if not args.cleanup:
            cmd.append("-no_cleanup")

        # Tell job control to copy compressed subjob's job dirs into the parent
        # job's dir when they complete:
        if args.SAVE:
            cmd.append("-SAVE")

        if args.fep_preparation:
            cmd.append("-fep_preparation")

        # NOTE: -HOST flag will be added automatically by JobDJ.

        debug("SUBMITTING CMD: %s" % cmd)
        jdj.addJob(cmd)

    if restarting:
        info("Restarting JobDJ instance...")
    else:
        info("Subjobs were prepared. Starting jobs...")

    total_steps = numjobs
    if backend:
        backend.setJobProgress(0, total_steps)

    # Just in case:
    jdj.dump(restart_file)

    # Dump the JobDJ instance to the restart file periodically:
    def callback():
        jdj.dump(restart_file)
        #info("Dumped JobDJ to file: %s" % restart_file)
        if backend:
            # Update the job progress after each completed job:
            steps = len(jdj.done_jobs)
            backend.setJobProgress(steps, total_steps)

    try:
        jdj.run(status_change_callback=callback)
    except RuntimeError as err:
        error(str(err))
        sys.exit(1)

    # Just in case:
    jdj.dump(restart_file)

    info("")  # Just a blank line

    # Check how many subjobs have completed and how many have failed:
    num_failed = len(jdj.failed_jobs)
    num_done = len(jdj.done_jobs)
    num_total = len(jdj.all_jobs)

    if num_failed == num_total:
        error("ERROR: All subjobs have failed; exiting the workflow.")
        sys.exit(1)

    info("%d of %d subjobs have completed successfully" % (num_done, num_total))
    if num_failed:
        msg = "WARNING: %d subjobs have failed! Their results will be " \
              "excluded from the merged results" % num_failed
        warning(msg)

    # Merge results from subjobs:
    combined_outfile, combined_refsfile = merge_subjob_results(
        main_jobname, sub_job_names, opt_affinity, fast_mode, args)

    if backend:
        # If running under job control
        backend.addOutputFile(combined_outfile)
        backend.addOutputFile(combined_refsfile)
        backend.setStructureOutputFile(combined_outfile)


def main(parser):
    """
    Run the main function using `parser` to set the options.
    Parses the arguments, and either runs the backend directly, or submits
    the backend under job control (startup script mode).

    :param parser: Parser used to get mutation options
    :type  parser: argparse.ArgumentParser

    """
    args = get_parser_args(parser)

    # This variable was set by the top-level script:
    JOBHOST = os.getenv('JOBHOST')  # first host (no ncpus)

    debug("JOBHOST: %s" % JOBHOST)

    use_jobcontrol = not args.NOJOBID
    debug("use_jobcontrol: %s" % use_jobcontrol)

    if use_jobcontrol:
        script_args = sys.argv[1:]
        submit_under_jobcontrol(args, script_args)

    else:
        # Run the backend (either job control was not being used, or we are
        # already running under job control.

        if args.NJOBS > 1:
            multijob_driver(args)
        else:
            # Either the user specified only 1 subjob OR we are already one
            # subjob
            run_backend(args)

        if (not args.subjob) and args.muts_file and (
                args.fast or args.mc and
            (args.refine_mut != 'prime_loop_prediction')):
            # If user submitted a muts_file, then reorder it per user's list.
            # Reordering should not be done on subjobs' output files
            # Reordering should not be run for jobs in discrete mode as these
            # already have ordered o/p and also have different format for
            # storing output cts.
            try:
                write_reordered_mutants(args.jobname + '-out.maegz',
                                        args.muts_file)
            except Exception as msg:
                error("ERROR: %s" % msg)
                sys.exit(1)


if __name__ == '__main__':

    parser = get_parser()
    main(parser)
