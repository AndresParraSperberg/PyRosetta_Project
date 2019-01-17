##################################################################
################### IMPORT MODULES ###############################
##################################################################

# import pyrosetta
from pyrosetta import *
init()

# import additional modules
from pyrosetta.teaching import *
import rosetta.protocols.rigid as rigid_moves
from rosetta.protocols.minimization_packing import *
from pyrosetta import PyMOLMover


# import and clean pdb (superantigen + TCR + MHC)
from pyrosetta.toolbox import cleanATOM
cleanATOM("2ICW_MHC.pdb")
pose = pose_from_pdb("2ICW_MHC.clean.pdb")

##################################################################
################### IMPORT POSE ##################################
##################################################################

# make full-atom starting pose
fa_starting = Pose()
fa_starting.assign(pose)
# make a full-atom working pose
fa_working = Pose()
fa_working.assign(pose)
# make a centroid pose
#switch = SwitchResidueTypeSetMover("centroid")
#switch.apply(pose)
#centroid = Pose()
#centroid.assign(pose)

##################################################################
######################### DESIGN #################################
##################################################################

# Strategy:
# 1. select residues on the surface

# pdb number D 55-72
print('this is what we want: ')
print(fa_working.pdb_info().pdb2pose('D', 55))
print(fa_working.pdb_info().pdb2pose('D', 72))
# desired range is 51 to 68 ()


# design only within 6A
#selectedResidues = pyrosetta.rosetta.core.select.residue_selector.ResidueRange()
#selectedResidues.set_start(51)
#selectedResidues.set_stop(68)

selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
selector.set_index('51-68')
print(selector.apply(pose))


#
#vec_test = pyrosetta.rosetta.core.select.get_residues_from_subset(selectedResidues)
#
#print(selectedResidues.start())
#print(selectedResidues.stop())
#print(vec_test)

sys.exit()




neighbors = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector() # not neighbors of those guys
neighbors.set_focus() # here give it your input residues as a boolean


neighbors = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(chain_A_cys_res, 6, True) # not neighbors of those guys
not_neighbors = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(neighbors) # a selector for all residues that are not chainA_cys_res


# Disable design
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
    pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), not_neighbors))

# Enable design on chain_A_cys_res
# THIS IS WHAT REPLACES THE RESFILE...
aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
aa_to_design.aas_to_keep("ADEFGHIKLMNPQRSTVWY")
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
    aa_to_design, chain_A_cys_res))

