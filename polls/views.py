# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from dynamic_rest import viewsets
from django.shortcuts import render
import os
from rest_framework import permissions
from rest_framework.response import Response
from rest_framework.permissions import AllowAny
from rest_framework.decorators import api_view, permission_classes, detail_route, list_route
from rest_framework import mixins
from rdkit import Chem
import random
from Bio.PDB.PDBParser import PDBParser
import time
from Bio.PDB.PDBIO import PDBIO
import residues_scanning_command
import prep_dock
import threading
from rosetta_workflow_all_scripts import rosetta_protein_prep, get_cst_file, change_pos, design_analysis
import rosetta_workflow_all_scripts
from django.utils.datastructures import MultiValueDictKeyError
import multiprocessing as mul
import time
from models import SubmitParamter, Onlinedock
from django.core.files import File
from email.mime.text import MIMEText
from email.header import Header
from smtplib import SMTP_SSL
from email.mime.multipart import MIMEMultipart
from email import encoders
from email.message import Message
from email.mime.base import MIMEBase
from dynamic_rest import viewsets
import serializers
from rest_framework.parsers import JSONParser
from polls.serializers import SubmitParamsSerializer, OnlinedockSerializer
from django.http import JsonResponse

import zipfile
import tempfile
def send_file_zipped(the_file, recipients, email_content, sender='1032847174@qq.com'):
    zf = tempfile.TemporaryFile(prefix='mail', suffix='zip')
    zip = zipfile.ZipFile(zf, 'w')
    zip.write(the_file)
    zip.close()
    zf.seek(0)

    ### Create the message
    themsg = MIMEMultipart()
    themsg['Subject'] = 'File %s' % the_file
    themsg['To'] = ', '.join(recipients)
    themsg['From'] = sender
    themsg.attach(MIMEText(email_content, 'html', 'utf-8'))
    themsg.preamble = 'I am not using a MIME-aware mail reader.\n'
    msg = MIMEBase('application', 'zip')
    msg.set_payload(zf.read())
    encoders.encode_base64(msg)
    msg.add_header('Content-Disposition', 'attachment', filename=the_file)
    themsg.attach(msg)
    themsg = themsg.as_string()

    ### Send the message
    import smtplib
    host_server = 'smtp.qq.com'
    sender_mail_addr = '1032847174@qq.com'
    pwd = 'utxfxpzcpsnzbbcc'

    smtp = SMTP_SSL(host_server)
    smtp.set_debuglevel(1)
    smtp.ehlo(host_server)
    smtp.login(sender_mail_addr, pwd)
    smtp.sendmail(sender, recipients, themsg)
    smtp.quit()

    # smtp = smtplib.SMTP()
    # smtp.connect()
    # smtp.sendmail(sender, recipients, themsg)
    # smtp.close()
def get_pov_value(file):
    f = open(file)
    lines = f.readlines()
    f.close()
    for line in lines:
        if line.startswith('1'):
            value = float(line.split('|')[1].strip())
            return value

def main(job_name, mutation_radius, pov_radius, pH, mutation_info_list, protein, ligand_name, ligand_resseq, chain_id, email_addr):
    """
    :param job_name:
    :param mutation_radius:
    :param pov_radius:
    :param pH:
    :param mutation_info_list:
    :param protein:
    :param ligand_name:
    :param ligand_resseq:
    :param chain_id:
    :return:
    """
    current_time = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
    print current_time
    log_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'log')
    wild_protein_name = protein.name
    # job_dir = os.path.join(log_dir, job_name + '_' + current_time)
    job_dir = os.path.join(log_dir, job_name)
    if not os.path.exists(job_dir):
        os.mkdir(job_dir)
    wild_protein_file = os.path.join(job_dir, wild_protein_name)
    protein_str = protein.read()
    prep_dock.save_to_file(wild_protein_file, protein_str)

    prepare_protein_name = wild_protein_name.split('.')[0] + '_prep.pdb'

    ### Prepare protein
    prep_dock.prep_protein(wild_protein_name, prepare_protein_name, job_dir, pH)

    ### make mutation
    prep_dock.get_mut_protein(job_dir, job_name, mut_info_list=mutation_info_list, mutation_radius=mutation_radius,
                    prepare_protein_name=prepare_protein_name)

    prepare_protein = os.path.join(job_dir, prepare_protein_name)
    mutation_protein_name = job_name + '_mut-2.pdb'
    mutation_protein = os.path.join(job_dir, mutation_protein_name)
    ### prep_pov
    prep_dock.get_pro_lig_povin((prepare_protein_name, prepare_protein, chain_id, ligand_resseq, ligand_name), pov_radius, protein_type='prep')
    prep_pov = os.path.join(job_dir, 'pov', 'prep', 'prep.log')
    ### mut_pov
    prep_dock.get_pro_lig_povin((mutation_protein_name, mutation_protein, chain_id, ligand_resseq, ligand_name), pov_radius, protein_type='mut')
    mut_pov = os.path.join(job_dir, 'pov', 'mut', 'mut.log')
    ### plip
    # prep_dock.get_plip_file(prepare_protein, mutation_protein)
    ### TMalign
    # prep_dock.TMalign(prepare_protein, mutation_protein)

    onlinedock, create = Onlinedock.objects.get_or_create(job_name=job_name)
    prep_protein_file = File(open(prepare_protein))
    mut_protein_file = File(open(mutation_protein))
    prep_pov_file = File(open(prep_pov))
    mut_pov_file = File(open(mut_pov))

    prep_pov_value = get_pov_value(prep_pov)
    mut_pov_value = get_pov_value(mut_pov)




    onlinedock.prep_protein.save(prepare_protein_name, prep_protein_file)
    onlinedock.mut_protein.save(mutation_protein_name, mut_protein_file)
    onlinedock.prep_pov.save('prep.log', prep_pov_file)
    onlinedock.mut_pov.save('mut.log', mut_pov_file)
    onlinedock.prep_pov_value = prep_pov_value
    onlinedock.mut_pov_value = mut_pov_value
    onlinedock.save()

    os.chdir(job_dir)
    os.system('zip related_info ' + prepare_protein + ' ' + mutation_protein + ' ' + prep_pov + ' ' + mut_pov)
    email_content = "Wellcome to Jianping Lin Group server~~"
    print(os.path.curdir)
    related_info = os.path.join(os.path.curdir, 'related_info.zip')
    send_file_zipped(related_info, email_addr, email_content=email_content)

def test():
    job_dir = '/home/jianping/django_test/longge/polls/log/1111/pov'


# job_id = 0
@api_view(['POST'])
@permission_classes([permissions.AllowAny])
def online_docking(request):
    # current_time = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
    # print current_time
    # log_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'log')

    job_name = request.data['job_name']
    mutation_radius = request.data['mutation_radius'] ### mutation radius
    pov_radius = str(request.data['pov_radius']) ### povelty radius
    pH = request.data['pH']
    mutation_info_list = request.data['mutation_info_list'] ### [chain, position, residue, ]
    protein = request.data['protein_file']
    ligand_name = request.data['ligand_name']
    ligand_resseq = int(request.data['ligand_resseq'])
    chain_id = request.data['chain_id']
    email_addr = request.data['email_addr']

    # main(job_name, mutation_radius, pov_radius, pH, mutation_info_list, protein, ligand_name, ligand_resseq, chain_id, email_addr)
    t = threading.Thread(target=main, args=(job_name, mutation_radius, pov_radius, pH, mutation_info_list, protein, ligand_name, ligand_resseq, chain_id, email_addr))
    t.setDaemon(False)
    t.start()
    return Response('Conguratulations, you have submitted successfully!!!')

@api_view(['POST'])
@permission_classes([permissions.AllowAny])
def prepare_protein(request):
    job_name = request.data['job_name']
    protein = request.data['protein']
    job_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Enzyme_design')
    work_dir = os.path.join(job_dir, job_name)
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)
    protein_name, protein_name_pure = protein.name, protein.name.split('.')[0]
    local_protein_file = os.path.join(work_dir, protein_name)
    protein_str = protein.read()
    prep_dock.save_to_file(local_protein_file, protein_str)

    os.chdir(work_dir)
    protein_renumber_name, protein_renumber = protein_name_pure + '_renumber', protein_name_pure + '_renumber.pdb'
    os.system(
        'python ../../rosetta_workflow_all_scripts/PDB_renumber.py -i ' + protein_name + ' -a -r > ' + protein_renumber_name + '.pdb')
    params, create = SubmitParamter.objects.get_or_create(job_name=job_name)
    prt = File(open(local_protein_file))
    prt_renumber = File(open(protein_renumber))

    params.protein_file.save(protein_name, prt)
    params.protein_renumber_file.save(protein_renumber, prt_renumber)
    params.save()
    # return Response(params)
    serializer = SubmitParamsSerializer(params)
    return JsonResponse(serializer.data, safe=False)
    # return Response('Successfully')

@api_view(['POST'])
@permission_classes([permissions.AllowAny])
def first_step(request):
    job_name = request.data['job_name']
    # protein = request.data['protein']
    ligand = request.data['ligand']
    other_ligands = request.data['other_ligands'] ### ['A','215','MG','A','218','HOH','A','217','ATP']
    other_ligands = other_ligands.split('[')[1].split(']')[0].split(',')
    other_ligands = [str(i) for i in other_ligands]
    res_chain = request.data['res_chain']  # 'A'
    res_ligand_chain = request.data['res_ligand_chain'] ## 'A'
    res_ligand_ID = request.data['res_ligand_ID'] ### '216'
    res_ligand_name = request.data['res_ligand_name'] ### 'PRP'
    # design_ligand_name = request.data['design_ligand_name'] ### 'ABC'

    ### third step ###
    # CST_A_chain_name = request.data['CST_A_chain_name']

    current_time = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
    print current_time
    job_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Enzyme_design')

    work_dir = os.path.join(job_dir, job_name)
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)

    # protein_name, protein_name_pure = protein.name, protein.name.split('.')[0]
    # local_protein_file = os.path.join(work_dir, protein_name)
    # protein_str = protein.read()
    # prep_dock.save_to_file(local_protein_file, protein_str)

    ligand_name, ligand_name_pure = ligand.name, ligand.name.split('.')[0]
    local_ligand_file = os.path.join(work_dir, ligand_name)
    ligand_str = ligand.read()
    prep_dock.save_to_file(local_ligand_file, ligand_str)
    os.chdir(work_dir)
    # protein_renumber_name, protein_renumber = protein_name_pure + '_renumber', protein_name_pure + '_renumber.pdb'
    # os.system('python ../../rosetta_workflow_all_scripts/PDB_renumber.py -i ' + protein_name + ' -a -r > ' + protein_renumber_name + '.pdb')
    os.system('python ../../rosetta_workflow_all_scripts/design_ligand_prep.py ' + ligand_name)
    while True:
        if os.path.exists(ligand_name_pure+'.params'):
            break
    os.system('cp ../../rosetta_workflow_all_scripts/match.flags ./')
    os.system('cp ../../rosetta_workflow_all_scripts/match_grid.flags ./')

    for filename in os.listdir(work_dir):
        if filename.endswith('renumber.pdb'):
            protein_renumber_name = filename.split('.pdb')[0]
            protein_renumber = filename
            break

    prep_pdb, prep_pdb_pure = protein_renumber_name + '_prep.pdb', protein_renumber_name + '_prep'

    rosetta_protein_prep.prep_protein(protein_renumber, prep_pdb, res_chain, './')
    rosetta_protein_prep.get_ligand(prep_pdb, res_ligand_chain, res_ligand_ID, res_ligand_name)

    ### my code ###
    step = 3
    other_ligands_class_list = [other_ligands[i: i+step] for i in range(0, len(other_ligands), step)]
    os.system('cp ' + protein_renumber_name + '_chain' + res_chain + '.pdb combi_ligands.pdb')

    if len(other_ligands) < 3:
        print 'There are no ligands that need to be retained'
        # os.system('cp ' + protein_renumber_name + '_chain' + res_chain + '.pdb combi_ligands.pdb')
    else:
        i = 0
        for cls in other_ligands_class_list:
            combi_name = '_'.join(cls)
            print combi_name
            rosetta_protein_prep.get_ligand(protein_renumber, cls[0], cls[1], cls[2])
            last_out_name = protein_renumber_name + '_chain' + combi_name + '.pdb'
            last_out_name_mol2 = protein_renumber_name + '_chain' + combi_name + '.mol2'
            rosetta_protein_prep.combi_pdb('combi_ligands.pdb', last_out_name)

            if cls[2] != 'HOH' and len(cls[2]) == 3:
                i += 1
                os.system('obabel -ipdb ' + last_out_name + ' -omol2 -O ' + last_out_name_mol2)
                os.system('python /home/jianping/Programs/rosetta_bin_linux_2018.09.60072_bundle/main/source/scripts/python/public/molfile_to_params.py ' + last_out_name_mol2 + '-n LG' + str(i))

    os.system("sed -i '/^TER/c'TER'' combi_ligands.pdb")
    rosetta_protein_prep.get_grid('../../rosetta_workflow_all_scripts/match_grid.flags', prep_pdb_pure, res_chain, res_ligand_chain, res_ligand_ID, res_ligand_name)
    rosetta_protein_prep.get_match_flags('../../rosetta_workflow_all_scripts/match.flags', res_chain, 'ABC', prep_pdb_pure, ligand_name_pure)
    os.system('/home/jianping/Programs/rosetta_bin_linux_2018.09.60072_bundle/main/source/bin/gen_lig_grids.linuxgccrelease -database /home/jianping/Programs/rosetta_bin_linux_2018.09.60072_bundle/main/database @match_grid_out.flags')
    os.system('cp ' + protein_renumber + ' ./renumber.pdb')


    ### update database ###
    # params, created = SubmitParamter.objects.get_or_create(
    #     job_name=job_name,
    #     other_ligands=other_ligands,
    #     res_chain=res_chain,
    #     res_ligand_chain=res_ligand_chain,
    #     res_ligand_name=res_ligand_name
    # )
    params = SubmitParamter.objects.get(job_name=job_name)
    params.other_ligands = other_ligands
    params.res_chain = res_chain
    params.res_ligand_chain = res_ligand_chain
    params.res_ligand_name = res_ligand_name

    # prt = File(open(local_protein_file))
    lgd = File(open(local_ligand_file))
    # prt_renumber = File(open(protein_renumber))
    ligand_params_file = File(open(ligand_name_pure+'.params'))
    pos_file = os.path.join('./inputs', prep_pdb_pure+'_chain'+res_chain+'.pdb_0.pos')
    pos_file_name = prep_pdb_pure+'_chain'+res_chain+'.pdb_0.pos'
    inputs_pos_file = File(open(pos_file))

    # params.protein_file.save(protein_name, prt)
    params.ligand_file.save(ligand_name, lgd)
    # params.protein_renumber_file.save(protein_renumber, prt_renumber)
    params.ligand_params_file.save(ligand_name_pure+'.params', ligand_params_file)
    params.inputs_pos_file.save(pos_file_name, inputs_pos_file)

    params.save()
    serializer = SubmitParamsSerializer(params)
    return JsonResponse(serializer.data, safe=False)
    # return Response('Successful')

@api_view(['POST'])
@permission_classes([permissions.AllowAny])
def second_step(request):
    job_name = request.data['job_name']
    constrain_info = request.data['constrain_info']  ### A:216:PRP:O2B:PB:O3A:A:131:ASP:OD2:CG:CB-O2B:PB:O3A-0.20:10.0:10.0:10.0:10.0:10.0-100.0:60.0:60.0:60.0:60.0:60.0-0:360.0:360.0:360.0:360.0:360.0-1:1:1:1:1:1, or A:216:PRP:O2B:PB:O3A:A:131:ASP:OD2:CG:CB-type:OH
    cat_ID = request.data['cat_ID']
    # cst1 = request.data['cst1']
    # cst2 = request.data['cst2']
    # cst3 = request.data['cst3']
    # three_atoms = request.data['three_atoms'] ### O2B:PB:O3A, type:OH
    # CST_A_chain_name = request.data['CST_A_chain_name'] ### 'A'
    # CST_A_residue_ID = int(request.data['CST_A_residue_ID']) ### '216'
    # CST_A_residue_name = request.data['CST_A_residue_name'] ### 'PRP'
    # Atom_A1 = request.data['Atom_A1'] ### 'O2B'
    # Atom_A2 = request.data['Atom_A2'] ### 'PB'
    # Atom_A3 = request.data['Atom_A3'] ### 'O3A'
    # CST_B_chain_name = request.data['CST_B_chain_name'] ### 'A'
    # CST_B_residue_ID = int(request.data['CST_B_residue_ID']) ### '131'
    # CST_B_residue_name = request.data['CST_B_residue_name'] ### 'ASP'
    # Atom_B1 = request.data['Atom_B1'] ### 'OD2'
    # Atom_B2 = request.data['Atom_B2'] ### 'CG'
    # Atom_B3 = request.data['Atom_B3'] ### 'CB'
    renumber_pdb = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Enzyme_design', job_name, 'renumber.pdb')
    work_dir = os.path.dirname(renumber_pdb)
    os.chdir(work_dir)

    ### my code ###
    #_______________________________________________________________
    constrain_info_list = [cst.split('-') for cst in constrain_info.split(',') if cst is not '']
    # for constrain_info in constrain_info_list:
    #     if len(constrain_info) == 2:

    parse = PDBParser(PERMISSIVE=1)
    structure = parse.get_structure('renumber.pdb', renumber_pdb)
    w = open('match.cst', 'w')
    w.write('# cst constraint descriptior for renumber.pdb' + '\n\n\n')
    w.write('# NOTE\n\n\n')

    for idx, cst_info in enumerate(constrain_info_list):
        cst_result = get_cst_file.measure_dist_angle_dihe_new(structure, idx, cst_info)
        w.writelines(cst_result)
    w.close()

    # get_cst_file.measure_dist_angle_dihe(renumber_pdb, 'renumber.pdb', constrain_info_list, 'match.cst')
    # ____________________________________________________________

    # get_cst_file.measure_dist_angle_dihe(renumber_pdb, 'renumber.pdb', [(CST_A_chain_name, CST_A_residue_ID, CST_A_residue_name,
    #                                                                     Atom_A1, Atom_A2, Atom_A3, CST_B_chain_name,
    #                                                                     CST_B_residue_ID, CST_B_residue_name, Atom_B1,
    #                                                                     Atom_B2, Atom_B3), ], 'match.cst')
    os.system('cp match.cst ./inputs')

    inputs_dir = os.path.join(work_dir, 'inputs')
    os.chdir(inputs_dir)
    for filename in os.listdir(inputs_dir):
        if filename.endswith('_0.pos'):
            pos = os.path.join(inputs_dir, filename)
            os.system('cp ' + pos + ' ./pos.bk')
            change_pos.change_pos(filename, cat_ID)
    params = SubmitParamter.objects.get(job_name=job_name)
    params.constrain_info = constrain_info
    params.cat_ID = cat_ID
    match_cst_file = File(open('match.cst'))
    params.match_cst_file.save('match.cst', match_cst_file)
    params.save()
    # for filename in os.listdir(inputs_dir):
    #     if filename.endswith('.params'):
    #         os.system('/home/jianping/Programs/rosetta_bin_linux_2018.09.60072_bundle/main/source/bin/CstfileToTheozymePDB.linuxgccrelease -database /home/jianping/Programs/rosetta_bin_linux_2018.09.60072_bundle/main/database -extra_res_fa ' + filename + ' -match:geometric_constraint_file match.cst')
    return Response('Successful')

@api_view(['POST'])
@permission_classes([permissions.AllowAny])
def third_step(request):
    job_name = request.data['job_name']
    # user_specialized_cst_file = request.data['user_specialized_cst_file']
    job_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Enzyme_design')
    work_dir = os.path.join(job_dir, job_name)
    os.chdir(work_dir)
    # if user_specialized_cst_file:
    #     cst_str = user_specialized_cst_file.read()
    #     user_defined_cst_file = os.path.join(work_dir, 'inputs', 'match.cst')
    #     prep_dock.save_to_file(user_defined_cst_file, cst_str)
    try:
        cst_file = request.data['cst_file']
        cst_str = cst_file.read()
        user_defined_cst_file = os.path.join(work_dir, 'inputs', 'match.cst')
        prep_dock.save_to_file(user_defined_cst_file, cst_str)
        params = SubmitParamter.objects.get(job_name=job_name)
        new_cst_file = File(open(user_defined_cst_file))
        params.user_defined_cst_file.save('match.cst', new_cst_file)
        params.save()
    except MultiValueDictKeyError:
        pass
    try:
        os.system('/home/jianping/Programs/rosetta_bin_linux_2018.09.60072_bundle/main/source/bin/match.linuxgccrelease @match_out.flags -database /home/jianping/Programs/rosetta_bin_linux_2018.09.60072_bundle/main/database')
        # params = SubmitParamter.objects.get(job_name=job_name)
        # UM_pdb_list = []
        # for filename in os.listdir(os.path.join(work_dir, 'inputs')):
        #     if filename.startswith('UM'):
        #         file = os.path.join(work_dir, 'inputs', filename)
        #         UM_pdb = File(open(file))
        UM_pdb_list = [filename for filename in os.listdir(os.path.join(work_dir, 'inputs')) if filename.startswith('UM')]
        params.UM_pdb_count = len(UM_pdb_list)
        params.save()
        # return Response('Successful, there are {} UM***.pdb'.format(len(UM_pdb_list)))
        serializer = SubmitParamsSerializer(params)
        return JsonResponse(serializer.data, safe=False)
    except:
        return Response('Failed, please check the constraint file and submit again !!!')

from functools import wraps
def timethis(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(func.__name__, end-start)
        return result
    return wrapper

@timethis
def design_comand(match_file):
    command = "/home/jianping/Programs/rosetta_bin_linux_2018.09.60072_bundle/main/source/bin/enzyme_design.linuxgccrelease @design_out.flags -database /home/jianping/Programs/rosetta_bin_linux_2018.09.60072_bundle/main/database -s "+ match_file + " -out:file:o " + match_file + "_DE.out > " + match_file + "_design.log"
    os.system(command)

# def get_design_params(ligand_name, params=('6.0', '8.0', '10.0', '12.0', '5')):
# #     """
# #     :param ligand_name:
# #     :param params: (6.0, 8.0, 10, 12.0, 5)
# #     :return:
# #     """
# #     command = ''
# #     os.system(command)

from functools import partial
def get_design_params(ligand_name, params=None): ### ligand_name not startswith('LG') endswith('params')
    if params is None:
        params = ('6.0', '8.0', '10.0', '12.0', '5')
        return partial(get_design_params, ligand_name)(params)
    # command = ''
    command = "sed -e 's/res_ligand_params_file/design_" + ligand_name + ".params/g' -e 's/enz_score.out/enz_score_" + ligand_name + ".out/g' -e 's/-cut1 6.0/-cut1 " + params[0] + "/g' -e 's/-cut2 10.0/-cut2 " + params[1] + "/g' -e 's/-cut3 15.0/-cut3 " + params[2] + "/g' -e 's/-cut4 20.0/-cut4 " + params[3] + "/g' -e 's/-nstruct 5/-nstruct " + params[4] + "/g' design.flags > design_out.flags"
    os.system(command)

def send_email(email_addr, email_content, result_file):
    host_server = 'smtp.qq.com'
    sender_mail_addr = '1032847174@qq.com'
    pwd = 'utxfxpzcpsnzbbcc'
    receiver_mail_addr = email_addr
    mail_content = email_content
    mail_title = "JianpingLin's email"
    msg = MIMEMultipart()
    msg['Subject'] = Header(mail_title, 'utf-8')
    msg['From'] = sender_mail_addr
    msg['To'] = Header('Receiver', 'utf-8')

    msg.attach(MIMEText(mail_content, 'html', 'utf-8'))
    # att1 = MIMEText(open(result_file).read(), 'base64', 'utf-8')
    att1 = MIMEText(open(result_file).read(), 'base64')
    # import zipfile
    # att1 = MIMEText(zipfile.ZipFile(result_file), 'base64', 'utf-8')
    att1['Content-Type'] = 'application/octet-stream'
    att1['Content-Disposition'] = 'attachment; filename="match_design.tar.gz"'
    msg.attach(att1)

    smtp = SMTP_SSL(host_server)
    smtp.set_debuglevel(1)
    smtp.ehlo(host_server)
    smtp.login(sender_mail_addr, pwd)
    smtp.sendmail(sender_mail_addr, receiver_mail_addr, msg.as_string())
    smtp.quit()


@api_view(['POST'])
@permission_classes([permissions.AllowAny])
def fourth_step(request):
    job_name = request.data['job_name']
    design_mini_range = request.data['design_mini_range']###
    user_email = request.data['user_email']

    # design_cst = request.data['design_cst']
    job_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Enzyme_design')
    work_dir = os.path.join(job_dir, job_name)
    os.chdir(work_dir)
    for filename in os.listdir(work_dir): ### ligand_name　必须是提交的ｍｏｌ２，不应该是ＬＧ.mol2
        if filename.endswith('params') and not filename.startswith('LG'):
            ligand_name = filename.split('.params')[0]
            break
    match_design_dir = os.path.join(work_dir, 'match_design')
    if not os.path.exists(match_design_dir):
        os.mkdir(match_design_dir)

    # if design_cst != '':
    #     cst_str = design_cst.read()
    #     user_design_cst_file = os.path.join(work_dir, 'match_design', 'design.cst')
    #     prep_dock.save_to_file(user_design_cst_file, cst_str)
    # else:
    #
    try:
        design_cst = request.data['design_cst']
        cst_str = design_cst.read()
        user_design_cst_file = os.path.join(work_dir, 'match_design', 'design.cst')
        prep_dock.save_to_file(user_design_cst_file, cst_str)
    except MultiValueDictKeyError:
        os.system('cp ./inputs/match.cst ./match_design/design.cst')
    finally:
        os.system('mv UM*match*.pdb ./match_design')
        os.system('cp ../../rosetta_workflow_all_scripts/design.flags ./')
        ###To DO###
        # command = "sed -e 's/res_ligand_params_file/design_" + ligand_name + ".params/g' -e 's/enz_score.out/enz_score_" + ligand_name + ".out/g' design.flags > design_out.flags"
        # get_design_params(ligand_name, tuple(design_mini_range.split(';')))
        ####TO DO###
        # os.system(command)

        if design_mini_range != '':
            #design_mini_range = req0uest.data['design_mini_range']
            tpl_mini_range = tuple(design_mini_range.split(';'))
            if len(tpl_mini_range) != 5:
                return Response('Please check that the "Designable Range, Repackable Range and Number of Outputs" exists.')
            else:
                get_design_params(ligand_name, tpl_mini_range)
        else:
            get_design_params(ligand_name)

        os.system("sed -r '/^PDB_ROTAMERS/d' " + ligand_name + ".params > match_design/design_" + ligand_name + ".params")
        os.system('cp design_out.flags ./match_design')
        match_dir = os.path.join(work_dir, 'match_design')
        os.chdir(match_dir)
        match_file_list = [filename for filename in os.listdir(match_dir) if filename.startswith('UM')]

        # design_comand(match_file_list[0])
        ###Post user###
        # pool = mul.Pool(5)
        # pool.map(design_comand, match_file_list)
        # pool.close()
        # pool.join()
    design_analysis.design_score(ligand_name, './')
    params = SubmitParamter.objects.get(job_name=job_name)
    params.user_email = user_email
    design_ligandname_out = 'design_' + ligand_name.split('.')[0] + '.out'

    file = File(open(design_ligandname_out))
    params.design_ligand_name_out.save(design_ligandname_out, file)
    params.save()

    # os.chdir(work_dir)
    # os.system('zip -r match_design.zip match_design')
    # os.system('tar czvf match_design.tar.gz UM*DE*.pdb')
    os.system('zip match_design UM*DE*.pdb ' + design_ligandname_out)

    email_content = "Welcome to Jianping Lin's group"
    match_design_file = os.path.join('./', 'match_design.zip')
    # send_email(email_addr=user_email, email_content=email_content, result_file=design_ligandname_out)
    # send_email(email_addr=user_email, email_content=email_content, result_file=match_design_file)
    # send_file_zipped(design_ligandname_out, ['1032847174@qq.com'])
    send_file_zipped(match_design_file, user_email, email_content=email_content)
    serializer = SubmitParamsSerializer(params)
    return JsonResponse(serializer.data, safe=False)
    # return Response('Successfully, this process needs ')

def get_analysis_params_dic(params):
    dic = {}
    temp_list = params.split(',')
    for param in temp_list:
        name, value = param.split(':')
        dic[name] = value
    return dic

@api_view(['POST'])
@permission_classes([permissions.AllowAny])
def fifth_step(request):
    job_name = request.data['job_name']
    analysis_params = request.data['analysis_params'] ### all_cst value < 0.9\nSR_2_interf_E_1_5:-9,
    # analysis_dict = get_analysis_params_dic(analysis_params) ### {all_cst:0.9,  SR_2_interf_E_1_5:-9}

    job_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Enzyme_design')
    work_dir = os.path.join(job_dir, job_name)
    os.chdir(work_dir)
    for filename in os.listdir(work_dir):
        if filename.endswith('params') and not filename.startswith('LG'):
            ligand_name = filename.split('.params')[0]
            break
    match_dir = os.path.join(work_dir, 'match_design')
    os.chdir(match_dir)
    design_analysis.design_filter(ligand_name, analysis_params.strip())
    # design_analysis.design_score(ligand_name, './')
    analysis_command = 'perl /home/jianping/Programs/rosetta_bin_linux_2018.09.60072_bundle/main/source/src/apps/public/enzdes/DesignSelect.pl -d ' + 'design_'+ligand_name+'.out' + ' -c ' + 'design_'+ligand_name+'.filter' + ' -tag_column last > filtered_designs_' + ligand_name +'.out'
    print analysis_command
    os.system(analysis_command)
    # serializer = SubmitParamsSerializer(params)
    # return JsonResponse(serializer.data, safe=False)
    return Response('Successfully')

class SubmitParamsViewSet(viewsets.DynamicModelViewSet):
    queryset = SubmitParamter.objects.all()
    serializer_class = serializers.SubmitParamsSerializer

class OnlinedockViewSet(viewsets.DynamicModelViewSet):
    queryset = Onlinedock.objects.all()
    serializer_class = serializers.OnlinedockSerializer
