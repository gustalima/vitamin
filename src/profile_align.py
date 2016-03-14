# -*- coding: utf-8 -*-
#Aligning sequences from PDB seg file and FASTA file

from modeller import *
import sys

pdbfile = sys.argv[1]
fastaseg = sys.argv[2]


log.verbose()
env = environ()
aln = alignment(env)
mdl = model(env, file=pdbfile)
aln.append_model(mdl, align_codes=pdbfile)  
aln.append(file=fastaseg+'.seg', align_codes=('target'))

# The as1.sim.mat similarity matrix is used by default:
aln.align(gap_penalties_1d=(-600, -400))
aln.write(file='profile.ali', alignment_format='PIR')
aln.write(file='profile.pap', alignment_format='PAP')


print "\n\nAlignment files created sucessfully!\n"

