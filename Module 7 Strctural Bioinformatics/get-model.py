# Step 4: model building
#
# This script should produce two models, 1fdx_my.B99990001.pdb and
# 1fdx_my.B99990002.pdb.

from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()
env = environ(rand_seed=-12312)  # To get different models from another script

# directories for input atom files
env.io.atom_files_directory = ['../atom_files']

a = automodel(env,
              alnfile='wnt1-6ry3.ali',      # alignment filename
              knowns=('6ry3'),    # codes of the templates
              sequence='wnt1',           # code of the target
              assess_methods=(assess.DOPE, assess.GA341))  #request GA341 assessment
a.starting_model= 1                 # index of the first model
a.ending_model  = 5             # index of the last model
                                    # (determines how many models to calculate)
a.deviation = 4.0                   # has to >0 if more than 1 model

a.make()                            # do homology modelling
