import pymol,sys
pymol.finish_launching()
from pymol import cmd
cmd.load('1vsn.pdb','1vsn')
#cmd.load('1vsn_2.pdb','1vsn2')

dst = cmd.get_distance('/1vsn//A/GLY`66/O','/1vsn//A/NFT`283/O23')
ang = cmd.get_angle('/1vsn//A/GLY`66/O','/1vsn//A/NFT`283/O23','/1vsn//A/NFT`283/C22')
dihe = cmd.get_dihedral('/1vsn//A/GLY`66/O','/1vsn//A/NFT`283/O23','/1vsn//A/NFT`283/C22','/1vsn//A/NFT`283/N26')
print "distance=""%.2f"%dst
print "angle=""%.2f"%ang
print "dihedral_angle=""%.2f"%dihe
#outfi.write("distance=%3f\n"%dst)
cmd.quit()
