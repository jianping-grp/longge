import os, sys
def get_conf(input_ligand,name,dir_path):
#	ligand_name = input_ligand.split('.')[0]
#	os.chdir(dir_path)
	os.system('mol2convert -imol2 '+name+'.mol2'+' -omae '+name+'.mae')
	print 'mol2convert -imol2 '+name+'.mol2'+' -omae '+name+'.mae'
	os.system('ligprep -g -s 1 -imae '+name+'.mae -omae '+name+'_ligprep.mae')
	print 'ligprep -g -s 1 -imae '+name+'.mae -omae '+name+'_ligprep.mae'
	while True:
		if os.path.exists(os.path.join(dir_path, name+'_ligprep.mae')):
			break
	os.system('confgenx -m 100 -no_cleanup '+name+'_ligprep.mae')
	while True:
		if os.path.exists(os.path.join(dir_path, name+'_ligprep-out.maegz')):
			break
	os.system('mol2convert -imae '+name+'_ligprep-out.maegz -omol2 '+name+'_ligprep-out_conf.mol2')
def get_params(res_ligand_name,ligand_name):
	# os.system('/home/liuc/programs/rosetta-3.5/rosetta_source/src/python/apps/public/molfile_to_params.py -n '+res_ligand_name+' -p '+ligand_name+' '+ligand_name+'_ligprep-out_conf.mol2')
	os.system('/home/chunfeng/Programs/rosetta_src_2018.33.60351_bundle/main/source/scripts/python/public/molfile_to_params.py -n '+res_ligand_name+' -p '+ligand_name+' '+ligand_name+'_ligprep-out_conf.mol2')

	os.system('cat '+ligand_name+'_????.pdb > '+ligand_name+'_confs.pdb')
	os.system('rm '+ligand_name+'_????.pdb')
	os.system('echo "PDB_ROTAMERS '+ligand_name+'_confs.pdb" >> '+ligand_name+'.params')


def main():
	input_ligand = sys.argv[1]
	ligand_path = './'
#	input_ligand = 'lig-high.mol2'
	ligand_name = input_ligand.split('.mol2')[0].split('/')[-1]
	print ligand_name
#	get_conf(input_ligand,'./')
	get_conf(input_ligand,ligand_name,ligand_path)
	res_ligand_name = 'ABC'
	get_params(res_ligand_name,ligand_name)

main()
