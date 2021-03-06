# -*- coding: utf-8 -*-
# Generated by Django 1.11 on 2018-11-08 11:48
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('polls', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='SubmitParamter',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('job_name', models.CharField(blank=True, max_length=1000, null=True)),
                ('protein_file', models.FileField(blank=True, null=True, upload_to='protein_file/')),
                ('ligand_file', models.FileField(blank=True, null=True, upload_to='ligand_file/')),
                ('other_ligands', models.CharField(blank=True, max_length=1000, null=True)),
                ('constrain_info', models.CharField(blank=True, max_length=20000, null=True)),
                ('cat_ID', models.CharField(blank=True, max_length=1000, null=True)),
                ('cst_file', models.FileField(blank=True, null=True, upload_to='cst_files/')),
                ('design_mini_range', models.CharField(blank=True, max_length=1000, null=True)),
                ('design_cst', models.FileField(blank=True, null=True, upload_to='design_cst_files/')),
                ('analysis_param', models.CharField(blank=True, max_length=3000, null=True)),
            ],
        ),
    ]
