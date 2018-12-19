import os

params_dict = {
    'all_cst': 'req all_cst value < ',
    'SR2': 'req SR_2_interf_E_1_5 value <'
}


def design_score(ligand_name, path):
	dirs = os.listdir(path)
	outfi=open('design_'+ligand_name+'.out','w')
	count = 0

	for fi in dirs:
		if fi.startswith('UM') and fi.split('_')[-1] == 'DE.out':
#			print fi
			count += 1
			infi = open(fi,'r')
			lines = infi.readlines()
			line_1 = lines[0]
			line_2 = lines[1:]
			if count == 1:
				for line in lines:
					outfi.write(line)
			if count > 1:
				for line in line_2:
					outfi.write(line)
			infi.close()
	outfi.close()
#	print outfi
	return outfi

def design_filter(ligand_name, param_str):
    des_fil = open('design_'+ligand_name+'.filter','w')

    param_list = param_str.split(';')
    des_fil.writelines(param_list)
    #des_fil.write(param_str)
    # for key, value in condition_dict.items():
    #     line = params_dict[key] + ' ' + value + '\n'
    #     des_fil.write(filter_condition)
    des_fil.close()



# def design_filter(ligand_name, condition_dict):
# 	des_fil = open('design_'+ligand_name+'.filter','w')
# 	for key, value in condition_dict.items():
#         line = params_dict[key] + ' ' + value + '\n'
#         des_fil.write(filter_condition)
#     des_fil.close()
#     return des_fil

if __name__ == '__main__':
    filter_condition = '''
    req all_cst value < 1.0
    req SR_2_interf_E_1_5 value < -10.0
    output sortmin total_score
    '''
    ligand_name = '56-73-5'
    path = './'
    #design_score(ligand_name, path)
    #design_filter(filter_condition)
    analysis_command = 'perl /home/liuc/programs/rosetta-3.5/rosetta_source/src/apps/public/enzdes/DesignSelect.pl -d ' + design_score(ligand_name, path).name + ' -c ' + design_filter(filter_condition).name + ' -tag_column last > filtered_designs_' + ligand_name +'.out'
    print analysis_command
    os.system(analysis_command)

