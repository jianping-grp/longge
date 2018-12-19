import os, sys
def prep_protein(input_pdb, prep_pdb, retained_chain, dir_path):
#	input_name = input_pdb.split('.pdb')[0].split('/')[-1]
#	prep_pdb = input_name+'_prep.pdb'
	os.system('prepwizard -watdist 0 -fillsidechains -noepik -s -propka_pH 7.5 -noimpref '+input_pdb+' '+ prep_pdb)
	#for example: prepwizard -watdist 0 -fillsidechains -noepik -s -propka_pH 7.5 -noimpref 1gm9.pdb 1gm9_prep.pdb
	while True:
		if os.path.exists(os.path.join(dir_path, prep_pdb)):
			os.system('sleep 5s')
			break
	os.system('python ../../rosetta_workflow_all_scripts/get_noHETATM_protein.py '+input_pdb+' '+retained_chain)
	os.system('python ../../rosetta_workflow_all_scripts/get_noHETATM_protein.py '+prep_pdb+' '+retained_chain)
	#for example: python2 ../get_noHETATM_protein.py 1gm9_prep.pdb AB

def get_ligand(prep_pdb,ligand_chain,ligand_id,ligand_name):
	os.system('python ../../rosetta_workflow_all_scripts/extract_ligand.py '+prep_pdb+' ' +ligand_chain+' ' +ligand_id+' '+ligand_name)
#	out_name = prep_pdb.split('.pdb')[0]+'_chain'+ligand_chain+'_'+ligand_id+'_'+ligand_name	
#	return out_name
	#for example: python22 extract_ligand.py 1gm9_prep.pdb B 1569 SOX
#def combi_pdb(retained_chain,ligand_chain,ligand_id,ligand_name):
def combi_pdb(first_pdb,second_pdb):
#	os.system('python22 change_chain_name_and_combi.py '+prep_pdb.split('.pdb')[0]+'_chain'+retained_chain+'.pdb'+' '+prep_pdb.split('.pdb')[0]+'_chain'+ligand_chain+'_'+ligand_id+'_'+ligand_name+'.pdb')
	os.system('python ../../rosetta_workflow_all_scripts/change_chain_name_and_combi.py '+first_pdb+' '+second_pdb)
	#for example: python2 ../change_chain_name_and_combi.py 1gm9_prep_chainAB.pdb 1gm9_prep_chainB_1568_CA.pdb
#	out_name = input_file1.split('.pdb')[0]+input_file2.split('prep')[-1]
#	return out_name
#	global combi_pdb_name = prep_pdb.split('.pdb')[0]+'_chain'+ligand_chain+'_'+ligand_id+'_'+ligand_name+'.pdb'

def get_grid(match_grid, prep_name, res_chain, res_ligand_chain, res_ligand_ID, res_ligand_name):
	if not os.path.exists('inputs'):
		os.mkdir('inputs')
	os.system('cp '+prep_name+'_chain'+res_ligand_chain+'_'+res_ligand_ID+'_'+res_ligand_name+'.pdb ./inputs')
	command = "sed -e 's/nohet.pdb/"+prep_name+"_chain"+res_chain+".pdb/g' -e 's/ligand.pdb/"+prep_name+"_chain"+res_ligand_chain+"_"+res_ligand_ID+"_"+res_ligand_name+".pdb/g' " + match_grid+" > ./match_grid_out.flags"
	os.system(command)
#	os.system('/home/liuc/programs/rosetta-3.5/rosetta_source/bin/gen_lig_grids.linuxgccrelease -database /home/liuc/programs/rosetta-3.5/rosetta_database @match_grid_out.flags')

def get_match_flags(match_flags, res_chain, design_ligand_name, prep_name, ligand_name_pure):
	if not os.path.exists('inputs'):
		os.mkdir('inputs')
	os.system('cp '+prep_name+'_chain'+res_chain+'.pdb ./inputs')
#	os.system('cp '+prep_name+'_chain'+res_chain+'.pdb_0.gridlig ./inputs')
#	os.system('cp '+prep_name+'_chain'+res_chain+'.pdb_0.pos ./inputs')
#	command = "sed -e 's/nohet.pdb/"+prep_name+"_chain"+res_chain+".pdb/g' -e 's/@@@/"+design_ligand_name+"/g' -e 's/grid_file/"+prep_name+"_chain"+res_chain+".pdb_0.gridlig/g' -e 's/list_pos_file/"+prep_name+"_chain"+res_chain+".pdb_0.pos/g' -e 's/res_ligand_params_file/"+ design_ligand_filename+".parmas/g' match.flags > match_out.flags"
	command = "sed -e 's/nohet.pdb/combi_ligands.pdb/g' -e 's/@@@/"+design_ligand_name+"/g' -e 's/grid_file/"+prep_name+"_chain"+res_chain+".pdb_0.gridlig/g' -e 's/list_pos_file/"+prep_name+"_chain"+res_chain+".pdb_0.pos/g' -e 's/res_ligand_params_file/" + ligand_name_pure +".params/g' match.flags > match_out.flags"

	os.system(command)
	os.system('cp match_out.flags ./inputs')
	os.system('cp combi_ligands.pdb ./inputs')
	os.system('cp match_grid_out.flags ./inputs')
	os.system('cp 138-08-9.params ./inputs')
	os.system('cp 138-08-9_confs.pdb ./inputs')
	#def cst2pdb()
#	os.system('/home/liuc/programs/rosetta-3.5/rosetta_source/bin/CstfileToTheozymePDB.linuxgccrelease @cst-to-pdb-flags -database /home/liuc/programs/rosetta-3.5/rosetta_database')


def main():
	prep_protein(input_pdb,res_chain,'./')
	get_ligand(prep_pdb,res_ligand_chain,res_ligand_ID,res_ligand_name)
	step = 0
	if len(other_ligands) < 3:
		print 'There are no ligands that need to be retained'
		os.system('cp '+input_name+'_chain'+res_chain+'.pdb combi_ligands.pdb')
	else:
		print 'Retained ligands : '+ str(other_ligands)
		for i in other_ligands:
			step = step +1
			if step <=3:
				if step % 3 == 0:
					combi_name = str(other_ligands[step-3:step]).split('[')[-1].split(']')[0]
					get_ligand(input_pdb,combi_name.split(',')[0],combi_name.split(',')[1].strip(),combi_name.split(',')[2].strip())
					last_out_name = input_pdb.split('.pdb')[0]+'_chain'+str(combi_name.split(',')[0]).replace('\'','')+'_'+str(combi_name.split(',')[1].strip()).replace('\'','')+'_'+str(combi_name.split(',')[2].strip()).replace('\'','')+'.pdb'
					combi_pdb(input_name+'_chain'+res_chain+'.pdb',last_out_name)

			if step > 3:
				if step % 3 == 0:
					combi_name = str(other_ligands[step-3:step]).split('[')[-1].split(']')[0]
					get_ligand(input_pdb,combi_name.split(',')[0],combi_name.split(',')[1].strip(),combi_name.split(',')[2].strip())
					last_out_name = input_pdb.split('.pdb')[0]+'_chain'+str(combi_name.split(',')[0]).replace('\'','')+'_'+str(combi_name.split(',')[1].strip()).replace('\'','')+'_'+str(combi_name.split(',')[2].strip()).replace('\'','')+'.pdb'
					combi_pdb('combi_ligands.pdb',last_out_name)
#	prep_protein(combi_pdb_name)
	os.system("sed -i '/^TER/c'TER'' combi_ligands.pdb")
#	get_conf(input_ligand,'./')
#	get_params()
	get_grid(match_grid)
	get_match_flags(match_flags)
	os.system('/home/jianping/Programs/rosetta_bin_linux_2018.09.60072_bundle/main/source/bin/gen_lig_grids.linuxgccrelease -database /home/jianping/Programs/rosetta_bin_linux_2018.09.60072_bundle/main/database @match_grid_out.flags')

if __name__ == '__main__':
    other_ligands = ['A','215','MG','A','218','HOH','A','219','HOH']
    match_grid = 'match_grid.flags'
    match_flags = 'match.flags'
    res_ligand_chain = 'A'
    res_ligand_ID = '216'
    res_ligand_name = 'PRP'
    res_chain = 'A'
    #input_ligand = 'lig-high.mol2'
    #ligand_name = input_ligand.split('.')[0]
    design_ligand_filename = '138-08-9'
    design_ligand_name = 'ABC'
    input_pdb = '1D6N_A68K_G102K_F35L_renumber.pdb'
    input_name = input_pdb.split('.pdb')[0].split('/')[-1]
    prep_pdb = input_name+'_prep.pdb'
    prep_name = prep_pdb.split('.pdb')[0]
    main()
    #os.system('/home/liuc/programs/rosetta-3.5/rosetta_source/bin/gen_lig_grids.linuxgccrelease -database /home/liuc/programs/rosetta-3.5/rosetta_database @match_grid_out.flags')