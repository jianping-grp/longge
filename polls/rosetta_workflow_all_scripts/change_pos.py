import os, sys

def change_pos(input_posfile, cat_ID):
    # command_bk = 'cp '+ input_posfile +' ' + input_posfile + '_bk'
    # os.system(command_bk)
    fi = open(input_posfile,'w')
    fi.write(cat_ID)
    fi.close()

if __name__ == '__main__':
    cat_ID = '131'
    input_posfile = sys.argv[1]
    change_pos(cat_ID)
