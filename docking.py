# Project Ideas:
# (1) Apply PyRosetta Docking to confirm native = lowest energy MHC binding structure
# (2) recapture binding orientation of native conformation (sampling)
# (3) See if protein has any regions that would be loaded onto MHC, if so remove and re-dock (in native conformation)
# note: for (3) we can just pick an arbitrary mutation, make it, and realign to native, dock_single and see energy

# import pyrosetta
from pyrosetta import *
init()

# import and clean pdb (superantigen + TCR + MHC)
from pyrosetta.toolbox import cleanATOM
cleanATOM("2ICW_MHC.pdb")
pose = pose_from_pdb("2ICW_MHC.clean.pdb")

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

# print pdb pose_tree
print(fa_working.fold_tree())

# change fold tree to keep DEF fixed and move H
from pyrosetta.teaching import *
setup_foldtree(fa_working, "DEF_H", Vector1([1])) # "DEF_H"
print(fa_working.fold_tree())

# see jump information
jump_num = 1
print(fa_working.jump(jump_num).get_rotation())
print(fa_working.jump(jump_num).get_translation())

# rotate and translate superantigen (8 degrees rot, 3 ang trans)
import rosetta.protocols.rigid as rigid_moves
pert_mover = rigid_moves.RigidBodyPerturbMover(jump_num, 8, 3) # .1 and 5 (interface object)
pert_mover.apply(fa_working)

# minimize the energy
scorefxn = get_fa_scorefxn()
from rosetta.protocols.minimization_packing import *
min_mover = MinMover()
movemap = MoveMap()
movemap.set_bb(True)
min_mover.movemap(movemap)
min_mover.score_function(scorefxn)
min_mover.apply(fa_working)

# see in PyMOL
from pyrosetta import PyMOLMover
pymol = PyMOLMover()
pymol.keep_history(True)
pymol.apply(fa_starting)
pymol.apply(fa_working)

# score pose
print('starting pose energy: ')
print(scorefxn(fa_starting))

print('working pose energy: ')
print(scorefxn(fa_working))


# testing
print(fa_working.pdb_info().pose2pdb(500))
