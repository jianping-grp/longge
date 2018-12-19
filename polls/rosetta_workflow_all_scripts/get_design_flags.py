import os, sys

def get_design_flags(ligand_name):
	
#	command = "sed -e 's/nohet.pdb/"+prep_name+"_chain"+res_chain+".pdb/g' -e 's/@@@/"+design_ligand_name+"/g' -e 's/grid_file/"+prep_name+"_chain"+res_chain+".pdb_0.gridlig/g' -e 's/list_pos_file/"+prep_name+"_chain"+res_chain+".pdb_0.pos/g' match.flags > match_out.flags"
	command = "sed -e 's/res_ligand_params_file/design_"+ligand_name+".params/g' -e 's/enz_score.out/enz_score_"+ligand_name+".out/g' design.flags > design_out.flags"
	os.system(command)
#	os.system('cp match_out.flags ./inputs')
	
#design_flags = 'design.flags'


##command
#home/liuc/programs/rosetta-3.5/rosetta_source/bin/enzyme_design.linuxgccrelease @design_out.flags -database /home/liuc/programs/rosetta-3.5/rosetta_database -s UM_15_D131_combi_ligands_match_1.pdb ##UM_***_***_combi_ligands_match_*.pdb
if __name__ == '__main__':
    ligand_name = "138-08-9"
    get_design_flags()
