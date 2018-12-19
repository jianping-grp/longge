# # -*- coding: utf-8 -*-
# from __future__ import unicode_literals
#
# from django.test import TestCase
#
# # Create your tests here.
from Bio.PDB.PDBParser import PDBParser
import os
base_dir = '/home/jianping/django_test/longge/polls/Enzyme_design/1111/'
pdb_file = os.path.join(base_dir, 'renumber.pdb')
parse = PDBParser(PERMISSIVE=1)
structure = parse.get_structure('renumber.pdb', pdb_file)


atom1 = structure[0]['A'][131]['OD2']
print atom1


