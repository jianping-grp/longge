import os

def split_pdb(pdb_name): #some .maegz can not change to .pdb automatic, if it can change automatic use if in this function ,else use else in this function
	if os.path.exists(pdb_name+'-2.pdb'):
		os.system('sed -i \'1d\' '+pdb_name+'-2.pdb')
		os.system('sed -i \'1i HEADER\t'+pdb_name+'\' '+pdb_name+'-2.pdb')  #add HEADER to pdb.
	else:
		infi = open(pdb_name+'.pdb')
		outfi = open(pdb_name+'-2.pdb','w')
		text = infi.read()
		infi.close()
		outfi.write(text.split('ENDMDL')[1]+'ENDMDL\nEND   ')
		outfi.close()
		os.system('sed -i \'1i HEADER\t'+pdb_name+'\' '+pdb_name+'-2.pdb')

def waiting(result_file):	#wait for the SCHRODINGER software completely
	while True:
		if os.path.exists(result_file+'-out.maegz'):
			os.system('pdbconvert -imae '+result_file+'-out.maegz -opdb '+result_file+'.pdb')  #change the schrodinger result maegz to pdb.
			break
		# else:
		# 	os.system('sleep 1m')


def main():

	os.system('$SCHRODINGER/run /home/jianping/Programs/schrodinger2015/mmshare-v30010/python/scripts/residue_scanning_backend.py -jobname 1d6n_mutate_chainA_ARG150ALA -residues A:150 -mut ALA -refine_mut prime_sidechain_bb -dist 5.0 1d6n_prep.pdb')
#check_file = '1d6n_mutate_chainA_ARG150ALA-out.maegz'
	waiting('1d6n_mutate_chainA_ARG150ALA')
	split_pdb('1d6n_mutate_chainA_ARG150ALA')

	os.system('$SCHRODINGER/run /home/jianping/Programs/schrodinger2015/mmshare-v30010/python/scripts/residue_scanning_backend.py -jobname 1d6n_mutate_chainAB_ARG150ALA -residues B:150 -mut ALA -refine_mut prime_sidechain_bb -dist 5.0 1d6n_mutate_chainA_ARG150ALA-2.pdb')
	waiting('1d6n_mutate_chainAB_ARG150ALA')
	split_pdb('1d6n_mutate_chainAB_ARG150ALA')

	os.system('$SCHRODINGER/run /home/jianping/Programs/schrodinger2015/mmshare-v30010/python/scripts/residue_scanning_backend.py -jobname 1d6n_mutate_chainAB_ARG150ALA_chainA_VAL157ALA -residues B:150 -mut ALA -refine_mut prime_sidechain_bb -dist 5.0 1d6n_mutate_chainAB_ARG150ALA-2.pdb')
	waiting('1d6n_mutate_chainAB_ARG150ALA_chainA_VAL157ALA')
	split_pdb('1d6n_mutate_chainAB_ARG150ALA_chainA_VAL157ALA')

	os.system('$SCHRODINGER/run /home/jianping/Programs/schrodinger2015/mmshare-v30010/python/scripts/residue_scanning_backend.py -jobname 1d6n_mutate_chainAB_ARG150ALA_chainAB_VAL157ALA -residues B:150 -mut ALA -refine_mut prime_sidechain_bb -dist 5.0 1d6n_mutate_chainAB_ARG150ALA_chainA_VAL157ALA-2.pdb')
	waiting('1d6n_mutate_chainAB_ARG150ALA_chainAB_VAL157ALA')
	split_pdb('1d6n_mutate_chainAB_ARG150ALA_chainAB_VAL157ALA')

if __name__ == '__main__':
    main()
