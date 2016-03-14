from modeller import *
import sys

alifile = sys.argv[1]
slct = sys.argv[2]

log.verbose()
env = environ()

env.libs.topology.read(file='$(LIB)/top_heav.lib')

# Read aligned structure(s):
aln = alignment(env)
aln.append(file='%s.ali' %alifile, align_codes='all')
aln_block = len(aln)

# Read aligned sequence(s):
aln.append(file='%s.seg' %slct, align_codes='target')

# Structure sensitive variable gap penalty sequence-sequence alignment:
aln.salign(output='', max_gap_length=20,
           gap_function=True,   # to use structure-dependent gap penalty
           alignment_type='PAIRWISE', align_block=aln_block,
           feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
           gap_penalties_1d=(-450, 0),
           gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
           similarity_flag=True)

aln.write(file='%s-mult.ali' %slct, alignment_format='PIR')
aln.write(file='%s-mult.pap' %slct, alignment_format='PAP')
