# import pyrosetta
from pyrosetta import *
init()

# import and clean pdb (superantigen + TCR + MHC)
from pyrosetta.toolbox import cleanATOM
cleanATOM("2ICW_MHC.pdb")
pose = pose_from_pdb("2ICW_MHC.clean.pdb")

############################################################

import rosetta.core.select.residue_selector as select

glycan_selector = select.GlycanResidueSelector()
residue_selector = select.ResidueIndexSelector()

my_selection = residue_selector.apply(pose)
