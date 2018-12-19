# coding=utf-8
import re

files = '/home/jianping/django_test/longge/polls/test_dir/1add_complex.pdb'
with open(files, 'r') as f:
    data = f.read().split('\n')
het = [n for n in data if n.startswith('HETATM')]


def get_chain(chain_lst, ligand_lst):
    reg = '|'.join(chain_lst)
    re_reg = re.compile(reg)
    reg_lig = '(%s)' % '|'.join(ligand_lst)
    re_reg_lig = re.compile(reg_lig)
    het = [n for n in data if n.startswith('HETATM')]
    ligand = [n for n in het if re_reg.match(n.split()[4]) and re_reg_lig.match(n.split()[3])]
    atom = [n for n in data if n.startswith('ATOM')]
    chains = [n for n in atom if re_reg.match(n.split()[4])]
    chains.extend(ligand)
    file_name = files.split('.')[0]+'_'+''.join(chain_lst)+'.pdb'
    with open(file_name, 'w') as w:
        w.write('\n'.join(chains))


def get_ligand(chain_lst, ligand_lst):
    reg = '(%s)' % '|'.join(ligand_lst)
    re_reg = re.compile(reg)
    reg_lig = '|'.join(chain_lst)
    re_reg_lig = re.compile(reg_lig)
    het = [n for n in data if n.startswith('HETATM')]
    ligand = [n for n in het if re_reg_lig.match(n.split()[4]) and re_reg.match(n.split()[3])]
    del_ligand = [n for n in het if n not in ligand]
    file_name_lig = files.split('.')[0]+'_'+'ligand'+'_and'.join(chain_lst)+'.pdb'
    with open(file_name_lig, 'w') as w:
        w.write('\n'.join(ligand))
    file_name_del_lig = files.split('.')[0] + '_chain' + ''.join(chain_lst) + '_del_' + '_and'.join(chain_lst) + '.pdb'
    with open(file_name_del_lig, 'w') as w:
        w.write('\n'.join(del_ligand))


def com_xyz(chain_lst, ligand_lst):
    het = [n for n in data if n.startswith('HETATM')]
    res = []
    for chain in chain_lst:
        reg = chain
        re_reg = re.compile(reg)
        for ligand in ligand_lst:
            reg_lig = ligand
            re_reg_lig = re.compile(reg_lig)
            ligand_ = [n for n in het if re_reg.match(n.split()[4]) and re_reg_lig.match(n.split()[3])]
            if len(ligand_) > 0:
                x, y, z = 0, 0, 0
                for n in ligand_:
                    x += float(n.split()[-6])
                    y += float(n.split()[-5])
                    z += float(n.split()[-4])
                x = x/len(ligand)
                y = y/len(ligand)
                z = z/len(ligand)
                res.append([chain, ligand, x, y, z])
    file_name = files.split('.')[0] + '.txt'
    with open(file_name, 'w') as w:
        for n in res:
            w.write('chain:%s, ligand:%s, 中心坐标为(%.3f,%.3f, %.3f)\n' % (n[0], n[1], n[2], n[3], n[4]))


get_chain(['A'], ['MOL'])
get_ligand(['A'], ['MOL'])
com_xyz(['A'], ['MOL'])

