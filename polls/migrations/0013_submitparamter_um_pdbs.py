# -*- coding: utf-8 -*-
# Generated by Django 1.11 on 2018-11-12 14:47
from __future__ import unicode_literals

import django.contrib.postgres.fields
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('polls', '0012_auto_20181109_1308'),
    ]

    operations = [
        migrations.AddField(
            model_name='submitparamter',
            name='UM_pdbs',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.FileField(blank=True, upload_to='UM_pdbs/'), blank=True, null=True, size=None),
        ),
    ]
