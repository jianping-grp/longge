from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import calc_angle, calc_dihedral
import math
from collections import namedtuple

def is_residue(residue):
    dic = ['ASP', 'ASN', 'GLY', 'ALA', 'ARG', 'ILE', 'CYS', 'GLU', 'GLN', 'HIS', 'LEU', 'MET', 'PHE', 'TRP', 'TRY', 'VAL', 'SER', 'PRO', 'LYS', 'THR']
    if residue in dic:
        return True
    else:
        return False
    # try:
    #     dic.get(residue)
    #     res_name = residue
    # except:
    #     res_name = 'H_'+residue
    # return res_name


# def measure_dist_angle_dihe(pdb_file, pdb_name, CST_list, output_file, cst1=(0.20, 10.0, 10.0,10.0,10.0,10.0), cst2=(100.0,60.0,60.0,60.0,60.0,60.0), cst3=(0,360.0,360.0,360.0,360.0,360.0), cst4=(1,1,1,1,1,1)):
def measure_dist_angle_dihe(pdb_file, pdb_name, CST_list, output_file):

    """
    CST_list = [(CST_A_chain_name, CST_A_residue_ID, CST_A_residue_name, Atom_A1,Atom_A2,Atom_A3,CST_B_chain_name,CST_B_residue_ID, CST_B_residue_name, Atom_B1,Atom_B2,Atom_B3
), ()]

    three_atom: ligand defined three atoms () or type:OH
    :param pdb_file:
    :param pdb_name:
    :param CST_A_residue_name:
    :param CST_A_chain_name:
    :param CST_A_residue_ID:
    :param Atom_A1:
    :param Atom_A2:
    :param Atom_A3:
    :param CST_B_chain_name:
    :param CST_B_residue_ID:
    :param CST_B_residue_name:
    :param Atom_B1:
    :param Atom_B2:
    :param Atom_B3:
    :param output_file:
    :return:
    """
    parse = PDBParser(PERMISSIVE=1)
    structure = parse.get_structure(pdb_name, pdb_file)
    w = open(output_file, 'w')
    w.write('# cst constraint descriptior for ' + pdb_name + '\n\n\n')
    w.write('# NOTE\n\n\n')
    CST = namedtuple('CST', ['CST_A_chain_name', 'CST_A_residue_ID', 'CST_A_residue_name', 'Atom_A1', 'Atom_A2', 'Atom_A3', 'CST_B_chain_name', 'CST_B_residue_ID', 'CST_B_residue_name', 'Atom_B1', 'Atom_B2', 'Atom_B3'])


    for idx, el in enumerate(CST_list, 1):
        cst, user_defined_atoms, cst1, cst2, cst3, cst4 = el
        cst, user_defined_atoms, cst1, cst2, cst3, cst4 = cst.split(':'), user_defined_atoms.split(':'), cst1.split(':'), cst2.split(':'), cst3.split(':'), cst4.split(':')
        cst = CST(*cst)
        atom_A1 = structure[0][cst.CST_A_chain_name][int(cst.CST_A_residue_ID)][cst.Atom_A1] if is_residue(cst.CST_A_residue_name) else structure[0][cst.CST_A_chain_name][('H_'+cst.CST_A_residue_name, int(cst.CST_A_residue_ID), ' ')][cst.Atom_A1]
        atom_A2 = structure[0][cst.CST_A_chain_name][int(cst.CST_A_residue_ID)][cst.Atom_A2] if is_residue(cst.CST_A_residue_name) else structure[0][cst.CST_A_chain_name][('H_'+cst.CST_A_residue_name, int(cst.CST_A_residue_ID), ' ')][cst.Atom_A2]
        atom_A3 = structure[0][cst.CST_A_chain_name][int(cst.CST_A_residue_ID)][cst.Atom_A3] if is_residue(cst.CST_A_residue_name) else structure[0][cst.CST_A_chain_name][('H_'+cst.CST_A_residue_name, int(cst.CST_A_residue_ID), ' ')][cst.Atom_A3]

        atom_B1 = structure[0][cst.CST_B_chain_name][int(cst.CST_B_residue_ID)][cst.Atom_B1] if is_residue(cst.CST_B_residue_name) else structure[0][cst.CST_B_chain_name][('H_'+cst.CST_B_residue_name, int(cst.CST_B_residue_ID), ' ')][cst.Atom_B1]
        atom_B2 = structure[0][cst.CST_B_chain_name][int(cst.CST_B_residue_ID)][cst.Atom_B2] if is_residue(cst.CST_B_residue_name) else structure[0][cst.CST_B_chain_name][('H_'+cst.CST_B_residue_name, int(cst.CST_B_residue_ID), ' ')][cst.Atom_B2]
        atom_B3 = structure[0][cst.CST_B_chain_name][int(cst.CST_B_residue_ID)][cst.Atom_B3] if is_residue(cst.CST_B_residue_name) else structure[0][cst.CST_B_chain_name][('H_'+cst.CST_B_residue_name, int(cst.CST_B_residue_ID), ' ')][cst.Atom_B3]

        # atom_A1 = structure[0][cst[0]][cst[1]][cst[3]] if is_residue(cst[2]) else structure[0][cst[0]][('H_'+cst[2], cst[1], ' ')][cst[3]]
        # atom_A2 = structure[0][cst[0]][cst[1]][cst[4]] if is_residue(cst[2]) else structure[0][cst[0]][('H_'+cst[2], cst[1], ' ')][cst[4]]
        # atom_A3 = structure[0][cst[0]][cst[1]][cst[5]] if is_residue(cst[2]) else structure[0][cst[0]][('H_'+cst[2], cst[1], ' ')][cst[5]]

        # print cst[6], cst[7], cst[8], cst[9]
        # print structure[0][cst[6]]
        # atom_B1 = structure[0][cst[6]][cst[7]][cst[9]] if is_residue(cst[8]) else structure[0][cst[6]][('H_'+cst[8], cst[7], ' ')][cst[9]]
        # atom_B2 = structure[0][cst[6]][cst[7]][cst[10]] if is_residue(cst[8]) else structure[0][cst[6]][('H_'+cst[8], cst[7], ' ')][cst[10]]
        # atom_B3 = structure[0][cst[6]][cst[7]][cst[11]] if is_residue(cst[8]) else structure[0][cst[6]][('H_'+cst[8], cst[7], ' ')][cst[11]]

        distanceAB = atom_A1 - atom_B1
        angleA = calc_angle(atom_A2.get_vector(), atom_A1.get_vector(), atom_B1.get_vector()) * 180 / math.pi### angle of A2, A1, B1
        angleB = calc_angle(atom_A1.get_vector(), atom_B1.get_vector(), atom_B2.get_vector()) * 180 / math.pi
        diheA = calc_dihedral(atom_A3.get_vector(), atom_A2.get_vector(), atom_A1.get_vector(), atom_B1.get_vector()) * 180 / math.pi
        diheAB = calc_dihedral(atom_A2.get_vector(), atom_A1.get_vector(), atom_B1.get_vector(), atom_B2.get_vector()) * 180 / math.pi
        diheB = calc_dihedral(atom_A1.get_vector(), atom_B1.get_vector(), atom_B2.get_vector(), atom_B3.get_vector()) * 180 / math.pi

        # w.write('# cst constraint descriptior for ' + pdb_name + '\n\n\n')
        # w.write('# NOTE\n\n\n')
        w.write('# block ' + str(idx) + ' for residue ABC' + ' ' + str(cst.CST_A_residue_name) + ' and residue ' + str(cst.CST_B_residue_ID) + ' ' + str(cst.CST_B_residue_name) + '\n\n')
        w.write('CST::BEGIN\n')

        # ligand_atom_type = user_defined_atoms.split(':')
        if 'type' not in user_defined_atoms:
            w.write('  TEMPLATE::   ATOM_MAP: 1 atom_name: '+' '.join(user_defined_atoms)+' '+'\n')
        else:
            w.write('  TEMPLATE::   ATOM_MAP: 1 atom_type: ' + ' '.join(user_defined_atoms[1:]) + ' \n')
        # w.write('  TEMPLATE::   ATOM_MAP: 1 residue3: ' + cst.CST_A_residue_name + '\n\n')
        w.write('  TEMPLATE::   ATOM_MAP: 1 residue3: ' + 'ABC' + '\n\n')
        w.write('  TEMPLATE::   ATOM_MAP: 2 atom_name: ' + cst.Atom_B1 + ' ' + cst.Atom_B2 + ' ' + cst.Atom_B3+' '+'\n')
        w.write('  TEMPLATE::   ATOM_MAP: 2 residue3: ' + cst.CST_B_residue_name + '\n\n')
        print cst1[0], cst2[0]
        w.write('  CONSTRAINT:: distanceAB: ' + "%6.2f" % float(distanceAB) + ' ' + "%4.2f" % float(cst1[0]) + ' ' + "%3.1f" % float(cst2[0]) + '   ' + str(int(cst3[0]))+ '    ' + str(int(cst4[0])) +'\n')
        w.write('  CONSTRAINT::    angle_A: ' + "%6.1f" % float(angleA) + ' ' + "%-4.1f" % float(cst1[1]) + ' ' + "%3.1f" % float(cst2[1]) + '  ' + "%-3.1f" % float(cst3[1]) + '  ' + str(int(cst4[1])) + '\n')
        w.write('  CONSTRAINT::    angle_B: ' + "%6.1f" % float(angleB) + ' ' + "%-4.1f" % float(cst1[2]) + ' ' + "%3.1f" % float(cst2[2]) + '  ' + "%-3.1f" % float(cst3[2]) + '  ' + str(int(cst4[2])) + '\n')
        w.write('  CONSTRAINT::  torsion_A: ' + "%6.1f" % float(diheA) + ' ' + "%-4.1f" % float(cst1[3]) + ' ' + "%3.1f" % float(cst2[3]) + '  ' + "%-3.1f" % float(cst3[3]) + '  ' + str(int(cst4[3])) + '\n')
        w.write('  CONSTRAINT:: torsion_AB: ' + "%6.1f" % float(diheAB) + ' ' + "%-4.1f" % float(cst1[4]) + ' ' + "%3.1f" % float(cst2[4]) + '  ' + "%-3.1f" % float(cst3[4]) + '  ' + str(int(cst4[4])) + '\n')
        w.write('  CONSTRAINT::  torsion_B: ' + "%6.1f" % float(diheB) + ' ' + "%-4.1f" % float(cst1[5]) + ' ' + "%3.1f" % float(cst2[5]) + '  ' + "%-3.1f" % float(cst3[5]) + '  ' + str(int(cst4[5])) + '\n')
        w.write('CST::END\n\n\n')
    w.close()

from functools import partial

# def measure_dist_angle_dihe(pdb_file, pdb_name, CST_list, output_file, cst1=(0.20, 10.0, 10.0,10.0,10.0,10.0), cst2=(100.0,60.0,60.0,60.0,60.0,60.0), cst3=(0,360.0,360.0,360.0,360.0,360.0), cst4=(1,1,1,1,1,1)):
def measure_dist_angle_dihe_new(structure, idx, CST_info):

    """
    CST_list = [(CST_A_chain_name, CST_A_residue_ID, CST_A_residue_name, Atom_A1,Atom_A2,Atom_A3,CST_B_chain_name,CST_B_residue_ID, CST_B_residue_name, Atom_B1,Atom_B2,Atom_B3
), ()]

    three_atom: ligand defined three atoms () or type:OH
    """
    result_list = []

    four_cst = ['0.20:10.0:10.0:10.0:10.0:10.0', '100.0:60.0:60.0:60.0:60.0:60.0', '0:360.0:360.0:360.0:360.0:360.0', '1:1:1:1:1:1']
    if len(CST_info) == 2:
        # return measure_dist_angle_dihe_new(structure, idx, CST_info.extend(four_cst))
        CST_info.extend(four_cst)
        return partial(measure_dist_angle_dihe_new, structure, idx)(CST_info)
    CST = namedtuple('CST', ['CST_A_chain_name', 'CST_A_residue_ID', 'CST_A_residue_name', 'Atom_A1', 'Atom_A2', 'Atom_A3', 'CST_B_chain_name', 'CST_B_residue_ID', 'CST_B_residue_name', 'Atom_B1', 'Atom_B2', 'Atom_B3'])
    cst, user_defined_atoms, cst1, cst2, cst3, cst4 = CST_info
    cst, user_defined_atoms, cst1, cst2, cst3, cst4 = cst.split(':'), user_defined_atoms.split(':'), cst1.split(
        ':'), cst2.split(':'), cst3.split(':'), cst4.split(':')
    cst = CST(*cst)
    atom_A1 = structure[0][cst.CST_A_chain_name][int(cst.CST_A_residue_ID)][cst.Atom_A1] if is_residue(
        cst.CST_A_residue_name) else \
    structure[0][cst.CST_A_chain_name][('H_' + cst.CST_A_residue_name, int(cst.CST_A_residue_ID), ' ')][cst.Atom_A1]
    atom_A2 = structure[0][cst.CST_A_chain_name][int(cst.CST_A_residue_ID)][cst.Atom_A2] if is_residue(
        cst.CST_A_residue_name) else \
    structure[0][cst.CST_A_chain_name][('H_' + cst.CST_A_residue_name, int(cst.CST_A_residue_ID), ' ')][cst.Atom_A2]
    atom_A3 = structure[0][cst.CST_A_chain_name][int(cst.CST_A_residue_ID)][cst.Atom_A3] if is_residue(
        cst.CST_A_residue_name) else \
    structure[0][cst.CST_A_chain_name][('H_' + cst.CST_A_residue_name, int(cst.CST_A_residue_ID), ' ')][cst.Atom_A3]

    atom_B1 = structure[0][cst.CST_B_chain_name][int(cst.CST_B_residue_ID)][cst.Atom_B1] if is_residue(
        cst.CST_B_residue_name) else \
    structure[0][cst.CST_B_chain_name][('H_' + cst.CST_B_residue_name, int(cst.CST_B_residue_ID), ' ')][cst.Atom_B1]
    atom_B2 = structure[0][cst.CST_B_chain_name][int(cst.CST_B_residue_ID)][cst.Atom_B2] if is_residue(
        cst.CST_B_residue_name) else \
    structure[0][cst.CST_B_chain_name][('H_' + cst.CST_B_residue_name, int(cst.CST_B_residue_ID), ' ')][cst.Atom_B2]
    atom_B3 = structure[0][cst.CST_B_chain_name][int(cst.CST_B_residue_ID)][cst.Atom_B3] if is_residue(
        cst.CST_B_residue_name) else \
    structure[0][cst.CST_B_chain_name][('H_' + cst.CST_B_residue_name, int(cst.CST_B_residue_ID), ' ')][cst.Atom_B3]

    distanceAB = atom_A1 - atom_B1
    angleA = calc_angle(atom_A2.get_vector(), atom_A1.get_vector(),
                        atom_B1.get_vector()) * 180 / math.pi  ### angle of A2, A1, B1
    angleB = calc_angle(atom_A1.get_vector(), atom_B1.get_vector(), atom_B2.get_vector()) * 180 / math.pi
    diheA = calc_dihedral(atom_A3.get_vector(), atom_A2.get_vector(), atom_A1.get_vector(),
                          atom_B1.get_vector()) * 180 / math.pi
    diheAB = calc_dihedral(atom_A2.get_vector(), atom_A1.get_vector(), atom_B1.get_vector(),
                           atom_B2.get_vector()) * 180 / math.pi
    diheB = calc_dihedral(atom_A1.get_vector(), atom_B1.get_vector(), atom_B2.get_vector(),
                          atom_B3.get_vector()) * 180 / math.pi
    result_list.append('# block ' + str(idx) + ' for residue ABC' + ' ' + str(cst.CST_A_residue_name) + ' and residue ' + str(cst.CST_B_residue_ID) + ' ' + str(cst.CST_B_residue_name) + '\n\n')
    result_list.append('CST::BEGIN\n')
    if 'type' not in user_defined_atoms:
        result_list.append('  TEMPLATE::   ATOM_MAP: 1 atom_name: '+' '.join(user_defined_atoms)+' '+'\n')
    else:
        result_list.append('  TEMPLATE::   ATOM_MAP: 1 atom_type: ' + ' '.join(user_defined_atoms[1:]) + ' \n')
    result_list.append('  TEMPLATE::   ATOM_MAP: 1 residue3: ' + 'ABC' + '\n\n')
    result_list.append('  TEMPLATE::   ATOM_MAP: 2 atom_name: ' + cst.Atom_B1 + ' ' + cst.Atom_B2 + ' ' + cst.Atom_B3+' '+'\n')
    result_list.append('  TEMPLATE::   ATOM_MAP: 2 residue3: ' + cst.CST_B_residue_name + '\n\n')
    result_list.append('  CONSTRAINT:: distanceAB: ' + "%6.2f" % float(distanceAB) + ' ' + "%4.2f" % float(cst1[0]) + ' ' + "%3.1f" % float(cst2[0]) + '   ' + str(int(cst3[0]))+ '    ' + str(int(cst4[0])) +'\n')
    result_list.append('  CONSTRAINT::    angle_A: ' + "%6.1f" % float(angleA) + ' ' + "%-4.1f" % float(cst1[1]) + ' ' + "%3.1f" % float(cst2[1]) + '  ' + "%-3.1f" % float(cst3[1]) + '  ' + str(int(cst4[1])) + '\n')
    result_list.append('  CONSTRAINT::    angle_B: ' + "%6.1f" % float(angleB) + ' ' + "%-4.1f" % float(cst1[2]) + ' ' + "%3.1f" % float(cst2[2]) + '  ' + "%-3.1f" % float(cst3[2]) + '  ' + str(int(cst4[2])) + '\n')
    result_list.append('  CONSTRAINT::  torsion_A: ' + "%6.1f" % float(diheA) + ' ' + "%-4.1f" % float(cst1[3]) + ' ' + "%3.1f" % float(cst2[3]) + '  ' + "%-3.1f" % float(cst3[3]) + '  ' + str(int(cst4[3])) + '\n')
    result_list.append('  CONSTRAINT:: torsion_AB: ' + "%6.1f" % float(diheAB) + ' ' + "%-4.1f" % float(cst1[4]) + ' ' + "%3.1f" % float(cst2[4]) + '  ' + "%-3.1f" % float(cst3[4]) + '  ' + str(int(cst4[4])) + '\n')
    result_list.append('  CONSTRAINT::  torsion_B: ' + "%6.1f" % float(diheB) + ' ' + "%-4.1f" % float(cst1[5]) + ' ' + "%3.1f" % float(cst2[5]) + '  ' + "%-3.1f" % float(cst3[5]) + '  ' + str(int(cst4[5])) + '\n')
    result_list.append('CST::END\n\n\n')
    return result_list




