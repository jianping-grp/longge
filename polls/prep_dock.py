from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO, Select
from collections import namedtuple
import os
import numpy as np
import command_TMalign_out
import residues_scanning_command

class ProteinSelect(Select):
    def __init__(self, chain_name, ligand_name, resseq):
        self.chain_name = chain_name
        self.ligand_name = ligand_name
        self.resseq = resseq

    def accept_residue(self, residue):
        if residue.get_parent().id == self.chain_name and residue.get_id() != ('H_'+self.ligand_name, self.resseq, ' '):
            return 1
        else:
            return 0

class LigandSelect(Select):
    def __init__(self, chain_name, ligand_name, resseq):
        self.chain_name = chain_name
        self.ligand_name = ligand_name
        self.resseq = resseq

    def accept_residue(self, residue):
        if residue.get_parent().id == self.chain_name and residue.get_id() == ('H_'+self.ligand_name, self.resseq, ' '):
            return 1
        else:
            return 0

def ligand_center(ligand):
    # print os.getcwd()
    ligand_dir, ligand_name = os.path.split(ligand)
    parser = PDBParser(PERMISSIVE=1)
    s = parser.get_structure(ligand_name, ligand)
    vector = np.zeros(3)
    for a in s.get_atoms():
        vector += a.get_coord()
        # print a.get_vector()
    return vector / len(list(s.get_atoms()))

# print(ligand_center('/home/jianping/django_test/longge/polls/A_PRO_only.pdb'))

def get_pov_in(ligand_center, pov_radius, job_dir, protein_name, protein_type):
    # print os.getcwd()
    base_dir = os.path.join(job_dir, 'pov')
    if not os.path.exists(base_dir):
        os.mkdir(base_dir)
    pov_dir = os.path.join(base_dir, protein_type)
    print ligand_center[1]
    ligand_center = '%.2f %.2f %.2f' % (ligand_center[0], ligand_center[1], ligand_center[2])
    if not os.path.exists(pov_dir):
        os.mkdir(pov_dir)
    os.chdir(pov_dir)
    os.system('cp ../../../../sample_POVME_input.ini .')
    command = "sed -e 's/xx.x yy.y zz.z R.R/" + ligand_center + " " + pov_radius + "/g' -e 's/INPUTFILE_NAME.PDB/\.\.\/\.\.\/" + protein_name + \
                   "/g' -e 's/.\/POVME_run\/POVME_9/\.\//g' sample_POVME_input.ini > POVME_output.ini"
    os.system(command)
    os.system('python ../../../../POVME2.py ./POVME_output.ini')
    # os.chdir(pov_dir)
    # prep_pov_dir = os.path.join(pov_dir, 'prep_pov')
    # mut_pov_dir = os.path.join(pov_dir, 'mut_pov')
    # if not os.path.exists(prep_pov_dir):
    #     os.mkdir(prep_pov_dir)
    # else:
    #     os.chdir(prep_pov_dir)
    #     os.system('cp ../../../../sample_POVME_input.ini .')
    #     prep_command = "sed -e 's/xx.x yy.y zz.z R.R/" + ligand_center + " " + pov_radius + "/g' -e 's/INPUTFILE_NAME.PDB/\.\.\/\.\.\/" + prep_protein_name + \
    #               "/g' -e 's/.\/POVME_run\/POVME_9/\.\//g' sample_POVME_input.ini > POVME_output.ini"
    #     os.system(prep_command)
    #     os.system('python ../../../../POVME2.py ./POVME_output.ini')
    #
    # if not os.path.exists(mut_pov_dir):
    #     os.mkdir(mut_pov_dir)
    # else:
    #     os.chdir(mut_pov_dir)
    #     os.system('cp ../../../../sample_POVME_input.ini .')
    #     mut_command = "sed -e 's/xx.x yy.y zz.z R.R/" + ligand_center + " " + pov_radius + "/g' -e 's/INPUTFILE_NAME.PDB/\.\.\/\.\.\/" + mut_protein_name + \
    #               "/g' -e 's/.\/POVME_run\/POVME_9/\.\//g' sample_POVME_input.ini > POVME_output.ini"
    #     os.system(mut_command)
    #     os.system('python ../../../../POVME2.py ./POVME_output.ini')


def get_plip_file(prepare_protein, mutate_protein):
    job_dir = os.path.split(prepare_protein)[0]
    plip_dir = os.path.join(job_dir, 'plip')
    if not os.path.exists(plip_dir):
        os.mkdir(plip_dir)
    os.system('cp ' + prepare_protein + ' ' + plip_dir)
    os.system('cp ' + mutate_protein + ' ' + plip_dir)
    os.chdir(plip_dir)
    os.system('plip -f ' + prepare_protein + '-tyxpo')
    os.system('plip -f ' + mutate_protein + '-tyxpo')

def TMalign(prepare_protein, mutate_protein):
    job_dir = os.path.split(prepare_protein)[0]
    TM_dir = os.path.join(job_dir, 'TMalign')
    if not os.path.exists(TM_dir):
        os.mkdir(TM_dir)
    os.chdir(job_dir)
    os.system('cp ../../../TMtools20170708/* ' + TM_dir)
    os.system('cp ' + prepare_protein + ' ' + TM_dir)
    os.system('cp ' + mutate_protein + ' ' + TM_dir)
    os.chdir(TM_dir)
    command_TMalign_out.TMalign(prepare_protein, mutate_protein, 'TM_out', './')

def get_pro_lig_povin(protein_info, pov_radius, protein_type):
    """
    :param protein_info: a tuple (protein_name, complex, chain_id, resseq, ligand_name)
    :return:
    """
    protein_name, complex, chain_id, resseq, ligand_name = protein_info
    job_dir = os.path.split(complex)[0]
    dock_protein_name = protein_type + '_chain_' + chain_id + '.pdb'
    dock_ligand_name = protein_type + '_ligand_' + chain_id + '_' + ligand_name + '_' + str(resseq) + '.pdb'
    dock_protein_dir = os.path.join(job_dir, dock_protein_name)
    dock_ligand_dir = os.path.join(job_dir, dock_ligand_name)
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(protein_name, complex)
    print(structure.get_chains())
    io = PDBIO()
    io.set_structure(structure)
    io.save(dock_protein_dir, ProteinSelect(str(chain_id), ligand_name, resseq))
    io.save(dock_ligand_dir, LigandSelect(str(chain_id), ligand_name, resseq))
    lig_center = ligand_center(dock_ligand_dir)
    print os.getcwd()
    # mut_protein_name = ''
    get_pov_in(lig_center, str(pov_radius), job_dir, dock_protein_name, protein_type)

# get_pro_lig_povin(('1d6n', '/home/jianping/django_test/longge/polls/log/1d6n_2018-09-10-13-57-17/1d6n_mut-2.pdb', 'A', 300, 'PPO'), 8.0, 'mut')

def prep_protein(wild_protein_name, prepare_protein_name, dir, pH=7.5):
    os.chdir(dir)
    os.system('prepwizard -watdist 0 -fillsidechains -s -propka_pH ' + pH + ' -r 0.3 -fix -f 2005 ' + wild_protein_name + ' ' + prepare_protein_name)
    while True:
        if os.path.exists(os.path.join(dir, prepare_protein_name)):
            break

def mutation_site(mut_info):
    mut_info_list = mut_info.split("[")[1].split(']')[0].split(',')
    step = 3
    result_list = [mut_info_list[i: i+step] for i in range(0, len(mut_info_list), step)]
    return result_list

# def get_mut_protein(job_dir, job_name=None, mut_info_list=None, mutation_radius=5.0, prepare_protein_name=None):
#     # mut_info_list = mut_info_list.split(')')[0].split('(')[1].split(', ')
#     mut_info_list = mutation_site(mut_info_list)
#     for idx, el in enumerate(mut_info_list):
#         chain = el[0]  ### wild type protein chain
#         position = el[1]  ### wild type protein position
#         residue = el[2]   ### mutate type protein residue name
#         os.chdir(job_dir)
#         os.system('$SCHRODINGER/run /home/jianping/Programs/Schrodinger_2018/mmshare-v4.1/python/scripts/residue_scanning_backend.py -jobname '
#                   + job_name + ' -residues ' + chain + ':' + str(position) + ' -mut ' + residue + ' -refine_mut prime_sidechain_bb -dist ' + str(mutation_radius) + ' ' + prepare_protein_name)
#         while True:
#             if os.path.exists(job_name+'-out.maegz'):
#                 os.system('pdbconvert -imae '+job_name+'-out.maegz -opdb ' + job_name+'_mut.pdb') ### better job_name+'_mut.pdb'
#                 break
#         residues_scanning_command.split_pdb(job_name)

def get_mut_protein(job_dir, job_name=None, mut_info_list=None, mutation_radius=5.0, prepare_protein_name=None):
    # mut_info_list = mut_info_list.split(')')[0].split('(')[1].split(', ')
    mut_info_list = mutation_site(mut_info_list)
    job_name = job_name + '_mut'
    for idx, el in enumerate(mut_info_list):
        # job_name = job_name + '_mut'
        chain = el[0]  ### wild type protein chain
        position = el[1]  ### wild type protein position
        residue = el[2]   ### mutate type protein residue name
        os.chdir(job_dir)
        if os.path.exists(job_name+'-2.pdb'):
            os.system(
                '$SCHRODINGER/run /home/jianping/Programs/Schrodinger_2018/mmshare-v4.1/python/scripts/residue_scanning_backend.py -jobname '
                + job_name + ' -residues ' + chain + ':' + str(
                    position) + ' -mut ' + residue + ' -refine_mut prime_sidechain_bb -dist ' + str(
                    mutation_radius) + ' ' + job_name+'-2.pdb')
            while True:
                if os.path.exists(job_name+'-out.maegz'):
                    os.system('pdbconvert -imae '+job_name+'-out.maegz -opdb ' + job_name + '.pdb') ### better job_name+'_mut.pdb'
                    os.remove(job_name+'-out.maegz')
                    break
            residues_scanning_command.split_pdb(job_name)
        else:
            os.system('$SCHRODINGER/run /home/jianping/Programs/Schrodinger_2018/mmshare-v4.1/python/scripts/residue_scanning_backend.py -jobname '
                  + job_name + ' -residues ' + chain + ':' + str(position) + ' -mut ' + residue + ' -refine_mut prime_sidechain_bb -dist ' + str(mutation_radius) + ' ' + prepare_protein_name)
            while True:
                if os.path.exists(job_name+'-out.maegz'):
                    os.system('pdbconvert -imae '+job_name+'-out.maegz -opdb ' + job_name + '.pdb') ### better job_name+'_mut.pdb'
                    os.remove(job_name+'-out.maegz')
                    break
            residues_scanning_command.split_pdb(job_name)
# get_mut_protein(job_dir='/home/jianping/django_test/longge/polls/log/1d6n_2018-09-10-13-57-17', job_name='1d6n', mut_info_list='[A, 150, ALA]', mutation_radius=5.0, prepare_protein_name='1d6n_prep.pdb')

def save_to_file(filename, contents):
    file = open(filename, 'w')
    file.write(contents)
    file.close()