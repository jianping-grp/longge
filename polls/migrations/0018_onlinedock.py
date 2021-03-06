# -*- coding: utf-8 -*-
# Generated by Django 1.11.16 on 2018-12-10 11:43
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('polls', '0017_submitparamter_um_pdb_count'),
    ]

    operations = [
        migrations.CreateModel(
            name='Onlinedock',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('job_name', models.CharField(blank=True, max_length=1000, null=True)),
                ('prep_protein', models.FileField(blank=True, null=True, upload_to=b'', verbose_name='origin_files/')),
                ('mut_protein', models.FileField(blank=True, null=True, upload_to=b'', verbose_name='mutate_files/')),
                ('prep_pov', models.FileField(blank=True, null=True, upload_to=b'', verbose_name='prep_pov_files/')),
                ('mut_pov', models.FileField(blank=True, null=True, upload_to=b'', verbose_name='mut_pov_files/')),
            ],
        ),
    ]
