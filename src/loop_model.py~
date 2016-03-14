# Loop refinement of an existing model
from modeller import *
from modeller.automodel import *
import sys

ini       = sys.argv[1]             #loop start
fin       = sys.argv[2]             #loop end    
bestmodel = sys.argv[3]+'.pdb'      #bestmodel pdb (from modeling routine)    
num_model = sys.argv[4]             #number of models to be build
seq_id    = sys.argv[5]             #id from sequence file

log.verbose()
env = environ()

# directories for input atom files
env.io.atom_files_directory = './:../atom_files'

# Create a new class based on 'loopmodel' so that we can redefine
# select_loop_atoms (necessary)
class MyLoop(loopmodel):
    # This routine picks the residues to be refined by loop modeling
    def select_loop_atoms(self):
        # 10 residue insertion 
        return selection(self.residue_range(ini, fin))

m = MyLoop(env,
           inimodel=bestmodel, # initial model of the target
           sequence=seq_id)          # code of the target

m.loop.starting_model= 1                # index of the first loop model 
m.loop.ending_model  = int(num_model)   # index of the last loop model
m.loop.md_level = refine.very_fast      # loop refinement method; this yields
                                        # models quickly but of low quality;
                                        # use refine.slow for better models

m.make()
