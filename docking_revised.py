##################################################################
########################### IMPORT ###############################
##################################################################

from pyrosetta import *
from pyrosetta.teaching import *
init("-mh:path:scores_BB_BB motif_dock/xh_16_ -mh:score:use_ss1 false -mh:score:use_ss2 false -mh:score:use_aa1 true -mh:score:use_aa2 true")

# load a protein:
pose = pose_from_pdb("2ICW_MHC.clean.pdb")
starting_pose = pose.clone()
working_pose = pose.clone()
native_pose = pose.clone()


##################################################################
################### SWITCH TO CENTROID ###########################
##################################################################

# switch to centroid
switch = SwitchResidueTypeSetMover("centroid")
switch.apply(working_pose)
switch.apply(native_pose)



##################################################################
################### CHANGE FOLD TREE #############################
##################################################################

# print pdb pose_tree
print(working_pose.fold_tree())

# change fold tree to keep DEF fixed and move H
setup_foldtree(working_pose, "DEF_H", Vector1([1])) # "DEF_H"
print(working_pose.fold_tree())

# set jump information
jump_num = 1
print(working_pose.jump(jump_num).get_rotation())
print(working_pose.jump(jump_num).get_translation())

##################################################################
#################### CENTROID ENERGY FUNCTION ####################
##################################################################

# this is the table used for centroid mode docking
motif_dock_score = create_score_function("motif_dock_score")

# what is the score?
print(motif_dock_score(working_pose)) # zero

##################################################################
############### BROAD LOCAL CENTROID MODE ########################
##################################################################

trans_mag = 3
rot_mag = 8
import rosetta.protocols.rigid as rigid_moves
rigid_body_mover = rigid_moves.RigidBodyPerturbMover(1, rot_mag, trans_mag)
rigid_body_mover.apply(working_pose)

kT = 0.8
mc = MonteCarlo(working_pose, motif_dock_score, kT)

import numpy as np
# note: change the magnitude dynamically inside...

# store the scores and Lrmsd
scores = np.array([])
Lrmsd = np.array([])

# current magnitudes
trans_mag = 3
rot_mag = 8
rigid_body_mover.trans_magnitude(trans_mag)
rigid_body_mover.rot_magnitude(rot_mag)

# add mover to the trial mover
trial_mover = TrialMover(rigid_body_mover, mc)
# trial_mover.apply(ras)

# loop through the monte carlo object
for i in range(10):
    for j in range(50):

        # apply the mover and update mc object
        trial_mover.apply(working_pose)
        #pymol.apply(working_pose)

        # store scores and Lrmsd
        scores = np.append(scores, motif_dock_score(working_pose))
        Lrmsd = np.append(Lrmsd, calc_Lrmsd(working_pose, native_pose, Vector1([1])))

    # get the acceptance rate
    acc = trial_mover.acceptance_rate()
    print(acc)


    # change the magnitude dynamically
    if acc < 0.5:
        trans_mag = trans_mag*0.9
        rot_mag = rot_mag*0.9
        rigid_body_mover.trans_magnitude(trans_mag)
        rigid_body_mover.rot_magnitude(rot_mag)

        # now, update the trial mover with the new rigid body mover
        trial_mover = TrialMover(rigid_body_mover, mc)

    elif acc > 0.5:
        trans_mag = trans_mag*1.1
        rot_mag = rot_mag*1.1
        rigid_body_mover.trans_magnitude(trans_mag)
        rigid_body_mover.rot_magnitude(rot_mag)

        # now, update the trial mover with the new rigid body mover
        trial_mover = TrialMover(rigid_body_mover, mc)

mc.show_scores()

##################################################################
##################### SWITCH BACK TO FULL-ATOM ###################
##################################################################

switch_back = SwitchResidueTypeSetMover("fa_standard")
switch_back.apply(working_pose)
#pymol.pymol_name("fa_pose")
#pymol.apply(working_pose)


##################################################################
############### FULL-ATOM ENERGY FUNCTION ########################
##################################################################

fa_sfxn = create_score_function("ref2015")


##################################################################
############### INTERFACE PACKER #################################
##################################################################

from pyrosetta.rosetta.protocols.interface import *
from pyrosetta.rosetta.core.select.residue_selector import *
from pyrosetta.rosetta.core.pack.task import *
from pyrosetta.rosetta.protocols.minimization_packing import *

def interface_packer(pose):
    """Create a packer to pack the interface

    Argument: pose whose interface needs to be packed

    return: a PackRotamersMover which repacks the interface,
    but doesn't change rotamers elsewhere
    """

    intf_resi_bool_vec = select_interface_residues(pose, "DEF_H", 8)

    intf_resi = ReturnResidueSubsetSelector(intf_resi_bool_vec)

    tf = TaskFactory()
    tf.push_back(operation.RestrictToRepacking()) # prevent design
    prevent_around_interface = operation.PreventRepackingRLT()
    # prevent repacking of all residues except interface residues
    # setting the 3rd argument to True flips the subset
    repack_only_intf_resi = operation.OperateOnResidueSubset(prevent_around_interface, intf_resi, True)
    tf.push_back(repack_only_intf_resi)

    packer = PackRotamersMover(fa_sfxn)
    packer.task_factory(tf)

    return packer


##################################################################
############### FULL-ATOM REFINEMENT #############################
##################################################################

# current magnitudes
trans_mag = 0.1
rot_mag = 5
rigid_body_mover.trans_magnitude(trans_mag)
rigid_body_mover.rot_magnitude(rot_mag)

# full-atom montecarlo object
kT = 0.8
mc_fa = MonteCarlo(working_pose, fa_sfxn, kT)

# trial mover
trial_mover = TrialMover(rigid_body_mover, mc_fa)


import numpy as np
# note: change the magnitude dynamically inside...

# store the scores and Lrmsd
scores = np.array([])
Lrmsd = np.array([])

# current magnitudes
trans_mag = 0.1
rot_mag = 5
rigid_body_mover.trans_magnitude(trans_mag)
rigid_body_mover.rot_magnitude(rot_mag)

# add mover to the trial mover
trial_mover = TrialMover(rigid_body_mover, mc_fa)
# trial_mover.apply(ras)

# loop through the monte carlo object
for i in range(5):
    for j in range(10):

        # apply the mover and update mc object
        trial_mover.apply(working_pose)
        #pymol.apply(working_pose)

        # store scores and Lrmsd
        scores = np.append(scores, motif_dock_score(working_pose))
        Lrmsd = np.append(Lrmsd, calc_Lrmsd(working_pose, native_pose, Vector1([1])))

    # get the acceptance rate
    acc = trial_mover.acceptance_rate()
    print(acc)

    # change the magnitude dynamically
    if acc < 0.5:
        trans_mag = trans_mag*0.9
        rot_mag = rot_mag*0.9
        rigid_body_mover.trans_magnitude(trans_mag)
        rigid_body_mover.rot_magnitude(rot_mag)

        # now, update the trial mover with the new rigid body mover
        trial_mover = TrialMover(rigid_body_mover, mc_fa)

    elif acc > 0.5:
        trans_mag = trans_mag*1.1
        rot_mag = rot_mag*1.1
        rigid_body_mover.trans_magnitude(trans_mag)
        rigid_body_mover.rot_magnitude(rot_mag)

        # now, update the trial mover with the new rigid body mover
        trial_mover = TrialMover(rigid_body_mover, mc_fa)

    # recover lowest observed pose
    mc_fa.recover_low(working_pose)

    # pack the rotamers at the interface
    interface_packer(working_pose)

mc_fa.show_scores()
mc_fa.recover_low(working_pose)
#pymol.apply(working_pose)


##################################################################
############### MINIMIZE #########################################
##################################################################

# Create a MinMover Object to minimize pose
from rosetta.protocols.minimization_packing import *
min_mover = MinMover()
mm = MoveMap()
mm.set_bb(True)
min_mover.movemap(mm)
min_mover.score_function(fa_sfxn)

min_mover.apply(working_pose)
