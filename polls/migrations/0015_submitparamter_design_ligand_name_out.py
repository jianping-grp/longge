# -*- coding: utf-8 -*-
# Generated by Django 1.11 on 2018-11-13 13:09
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('polls', '0014_remove_submitparamter_um_pdbs'),
    ]

    operations = [
        migrations.AddField(
            model_name='submitparamter',
            name='design_ligand_name_out',
            field=models.FileField(blank=True, null=True, upload_to='design_ligandname_out_files/'),
        ),
    ]
