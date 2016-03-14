from modeller import *
from modeller.scripts import complete_pdb
import sys


slctmod = sys.argv[1] #bestmodel
code    = sys.argv[2] #PDB template


log.verbose()    # request verbose output
env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

# read model file
mdl = complete_pdb(env, slctmod+'.pdb')

# Assess with DOPE:
s = selection(mdl)   # all atom selection
s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='%s.profile'%slctmod,
	      normalize_profile=True, smoothing_window=15)



# read model file


mdlt = complete_pdb(env, '%s.pdb' % code)

# Assess with DOPE:
s = selection(mdlt)   # all atom selection
s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='%s.profile' % code,
	      normalize_profile=True, smoothing_window=15)
