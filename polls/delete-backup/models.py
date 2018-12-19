# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models
from django.contrib.postgres.fields import ArrayField

class SubmitParamter(models.Model):

    ### first_step
    job_name = models.CharField(max_length=1000, blank=True, null=True)
    protein_file = models.FileField(upload_to='protein_files/', blank=True, null=True)
    protein_renumber_file = models.FileField(upload_to='protein_files/', blank=True, null=True)
    ligand_file = models.FileField(upload_to='ligand_files/', blank=True, null=True)
    other_ligands = models.CharField(max_length=1000, blank=True, null=True)
    res_chain = models.CharField(max_length=20, blank=True, null=True)
    res_ligand_chain = models.CharField(max_length=50, blank=True, null=True)
    res_ligand_name = models.CharField(max_length=50, blank=True, null=True)

    ligand_params_file = models.FileField(upload_to='ligand_params_files/', blank=True, null=True)
    inputs_pos_file = models.FileField(upload_to='inputs_pos_files/', blank=True, null=True)
    ###

    ### second_step ###
    constrain_info = models.CharField(max_length=20000, blank=True, null=True)
    cat_ID = models.CharField(max_length=1000, blank=True, null=True)
    match_cst_file = models.FileField(upload_to='match_cst_files/', blank=True, null=True)
    ###

    ### third_step ###
    user_defined_cst_file = models.FileField(upload_to='user_defined_cst_files/', blank=True, null=True)
    UM_pdb_count = models.IntegerField(null=True, blank=True)
    # UM_pdbs = ArrayField(models.FileField(upload_to='UM_pdbs/', blank=True), blank=True, null=True)
    ###

    ### fourth step ###
    design_mini_range = models.CharField(max_length=1000, blank=True, null=True)
    design_cst = models.FileField(upload_to='design_cst_files/', blank=True, null=True)
    design_ligand_name_out = models.FileField(upload_to='design_ligandname_out_files/', blank=True, null=True)
    user_email = models.EmailField(blank=True, null=True)
    ###

    ### fifth step ###
    analysis_param = models.CharField(max_length=3000, blank=True, null=True)
