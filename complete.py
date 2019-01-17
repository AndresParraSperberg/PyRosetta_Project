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
switch = SwitchResidueTypeSetMover("centroid")
switch.apply(pose)
centroid = Pose()
centroid.assign(pose)

##################################################################
################### RELAX PDB ####################################
##################################################################




##################################################################
################### MUTATE PROTEIN ###############################
##################################################################

# make score function
scorefxn = get_fa_scorefxn()

# mutate protein
print(fa_working.pdb_info().pdb2pose('H', 12))
mutater = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()
mutater.set_target(386)
mutater.set_res_name('PRO')
mutater.apply(fa_working)

# see in PyMOL
pymol = PyMOLMover()
pymol.pymol_name('mutated_superantigen')
pymol.apply(fa_working)


##################################################################
################### RELAX LOCALLY ################################
##################################################################

# relax pose, store lowest energy
# (i.e. import fast relax method from lecture7)
# task factory, what to design, where and how...
# (dock mcm protocol, get interface residues)
# interface energy optimization

#from fast_relax_function import *
#fast_relax(fa_working, fast_relax_rounds=1)

##########################
# CreateGlycanSequonMover
# Glycosylate mover
# man5
# glycan tree modeler
# sugarcoat
##########################

##################################################################
################### DOCK PROTEIN #################################
##################################################################

# initialize minimization
min_mover = MinMover()
movemap = MoveMap()
movemap.set_bb(True)
min_mover.movemap(movemap)
min_mover.score_function(scorefxn)

# print pdb pose_tree
print(fa_working.fold_tree())

# change fold tree to keep DEF fixed and move H
setup_foldtree(fa_working, "DEF_H", Vector1([1])) # "DEF_H"
print(fa_working.fold_tree())

# set jump information
jump_num = 1
print(fa_working.jump(jump_num).get_rotation())
print(fa_working.jump(jump_num).get_translation())

# display starting energy
starting_score = scorefxn(fa_starting)
print('starting pose energy: ')
print(starting_score)

# Loop through and store lowest energy docking pose
lowest_energy_pose = Pose()
lowest_energy = float('inf')

for i in range(3):
  # rotate and translate superantigen (8 degrees rot, 3 ang trans)
  pert_mover = rigid_moves.RigidBodyPerturbMover(jump_num, 8, 3)
  pert_mover.apply(fa_working)

  # minimize the energy
  min_mover.apply(fa_working)

  # score pose
  print('working pose energy: ')
  curr_score = scorefxn(fa_working)
  print(curr_score)

  # store values
  if curr_score < lowest_energy:
    lowest_energy = curr_score
    lowest_energy_pose.assign(fa_working)

##################################################################
################### SCORE ########################################
##################################################################

# find score of best protein
print('best energy: ')
print(scorefxn(lowest_energy_pose))

