import os

def TMalign(pdb1,pdb2,outfi,path_dir):
	os.system('./TMalign '+pdb1+' '+pdb2+' -o align_'+outfi)
	path = path_dir
	for file in os.listdir(path):
		if file.startswith('align_'+outfi):
			newname = file+'.pdb'
			os.rename(os.path.join(path,file),os.path.join(path,newname))
		else:
			continue

def main():
    TMalign('../1bzy_prep.pdb','../1d6n_prep.pdb','test1','./')

if __name__ == '__main__':
    main()
