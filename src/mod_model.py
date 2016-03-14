from modeller import *
from modeller.automodel import *

import sys

number = sys.argv[1]
alifile = sys.argv[2]



log.verbose()
env = environ()
templ = alifile.split('-')[0]

a = automodel(env, alnfile=alifile,
              knowns=templ, sequence='target',
              assess_methods=(assess.DOPE, assess.GA341))
a.starting_model= 1                 # index of the first model 

a.ending_model  = int(number)       # index of the last model

a.final_malign3d = 1                              
a.deviation = 4.0			    
a.max_var_iterations = 500 	    # maximal numb of iterations for the cycles
                                    # of the variable target function method
a.md_level = refine.slow    	    # what kind of optimization is done after
                              	    # the variable target function method:
                              	    # 'none'                ... nothing;
	                            # 'refine.very_fast'    ... very fast MD annealing;
	                            # 'refine.fast'         ... fast MD annealing;
	                            # 'refine.slow'           
			            # 'refine.very_slow'    ... thorough MD annealing;
			            
a.repeat_optimization = 4 	    #	how many times the whole optimization
	                            # schedule (variable target function
	                            # method and refinement) is repeated
	                            # for each initial model;

a.max_molpdf = 100E3		    # abort optimization of the current model if 
                              	    # the molecular pdf is larger than this and
	                            # continue with the next model;		    

a.write_intermediates = 0 	    # 0 ... do not write out intermediate atom files during optimization;
                                    # 1 ... write out intermediate atom files;	

a.final_malign3d = 1     	    # 0 ... do not do MALIGN3D and write 
                                    #       superposed templates & models
                                    #       at the end of 'model'
                                    # 1 ... do that


a.generate_method = generate.transfer_xyz  # how to build the initial model: 
                                    # 'generate_xyz' from internal coordinates 
                                    #                and write them to a file;
                                    # 'transfer_xyz' from template coordinates
                                    #                and write them to a file;
                                    # 'read_xyz'     read coordinates from 
                                    #                a file;

a.rand_method = randomize.xyz     # a method to perturb the initial model:
                                    # 'randomize_dihedrals' ... uses DEVIATION
                                    #                           in degrees;
                                    # 'randomize_xyz'       ... uses DEVIATION
                                    #                           in angstroms;
                                    # 'none'
			            
a.make()


        
