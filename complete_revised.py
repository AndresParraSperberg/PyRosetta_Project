##################################################################
##################################################################
##################################################################
##################################################################
############################# PROJECT ############################
##################################################################
##################################################################
##################################################################
##################################################################
# - R. Andres Parra Sperberg, Alican G, Charlotte



##################################################################
################### IMPORT MODULES ###############################
##################################################################

# import pyrosetta
from pyrosetta import *
init("-mh:path:scores_BB_BB motif_dock/xh_16_ -mh:score:use_ss1 false -mh:score:use_ss2 false -mh:score:use_aa1 true -mh:score:use_aa2 true") # for docking

# import additional modules
from pyrosetta.teaching import *
import rosetta.protocols.rigid as rigid_moves
from rosetta.protocols.minimization_packing import *
from pyrosetta import PyMOLMover
import numpy as np
from rosetta.protocols import minimization_packing as pack_min
from rosetta.protocols.relax import FastRelax

# import and clean pdb (superantigen + TCR + MHC)
from pyrosetta.toolbox import cleanATOM
cleanATOM("2ICW_MHC.pdb")
pose = pose_from_pdb("2ICW_MHC.clean.pdb")

##################################################################
############### FUNCTIONS ########################################
##################################################################

# 1. Interface Packer
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
################### IMPORT POSE ##################################
##################################################################

# make full-atom starting pose
fa_starting = Pose()
fa_starting.assign(pose)
# make a full-atom working pose
fa_working = Pose()
fa_working.assign(pose)

##################################################################
################### CONSTRAINT RELAX #############################
##################################################################

# Set up task Factory
tf_relax = TaskFactory()
tf_relax.push_back(operation.RestrictToRepacking())

# Set up packer for pose
packer_relax = pack_min.PackRotamersMover()
packer_relax.task_factory(tf_relax) # sets the task factory to the packer
packer_relax.apply(fa_working) # applies packer to pose

# Initializes fast relax class
fr = FastRelax()

# score function for fast relax
sfxn_relax = get_score_function()
fr.set_scorefxn(sfxn_relax)

# apply fast relax to pose
fr.apply(fa_working)

##################################################################
##################################################################
##################################################################
##################################################################
############################# DESIGN #############################
##################################################################
##################################################################
##################################################################
##################################################################



##################################################################
################### RESIDUE SELECTOR #############################
##################################################################

#find the pose numbering from pdb numbering
start_residue = fa_working.pdb_info().pdb2pose('D', 55)
end_residue = fa_working.pdb_info().pdb2pose('D', 72)

#first define the main helix of the antigen in contact with the protein
antigen_helix_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
antigen_helix_selector.set_index("{}-{}".format(start_residue, end_residue))
antigen_helix_residue_vector = antigen_helix_selector.apply(fa_working)

#find the neighbors
antigen_helix_neighbor_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(antigen_helix_residue_vector, 8, True)
antigen_helix_neighbors = antigen_helix_neighbor_selector.apply(fa_working)
not_antigen_helix_neighbors = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(antigen_helix_neighbor_selector)


#define the chain H residues belonging to the protein
H_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('H')
DF_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('D,E,F')

#combine the chain selectors with the region selectors
combine_H_selectors = pyrosetta.rosetta.core.select.residue_selector.AND_combine(antigen_helix_neighbor_selector, H_selector)
not_combine_H_selectors = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(combine_H_selectors)
H_residues_to_design = combine_H_selectors.apply(fa_working)

combine_DF_selectors = pyrosetta.rosetta.core.select.residue_selector.AND_combine(antigen_helix_neighbor_selector, DF_selector)
not_combine_DF_selectors = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(combine_DF_selectors)
DF_residues_to_repack = combine_DF_selectors.apply(fa_working)


##################################################################
########################## PACKER ################################
##################################################################

# setup the packer-task
tf_design = pyrosetta.rosetta.core.pack.task.TaskFactory()

# These are pretty standard
tf_design.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
tf_design.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
tf_design.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())


# Disable repacking/design on non-interface
tf_design.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
    pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(), not_antigen_helix_neighbors))

# Disable design (just repack) on MHC
tf_design.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
    pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), DF_residues_to_repack))

# Convert the task factory into a PackerTask
packer_task_design = tf_design.create_task_and_apply_taskoperations(fa_working)

# View the PackerTask
print(packer_task_design)


##################################################################
########################## MINIMIZER #############################
##################################################################

# setup the minimizer
mm_design = pyrosetta.rosetta.core.kinematics.MoveMap()
mm_design.set_bb(True)
mm_design.set_chi(True)
mm_design.set_jump(True) # psuedo docking...


##################################################################
######################## FAST DESIGN #############################
##################################################################

# setup fast design
scorefxn_design = pyrosetta.create_score_function("ref2015_cart.wts")
fast_design = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign(scorefxn_in=scorefxn_design, standard_repeats=1)
fast_design.cartesian(True)
fast_design.set_task_factory(tf_design)
fast_design.set_movemap(mm_design)
fast_design.set_relaxscript(1) # lower it to 2....
#fast_design.minimize_bond_angles(True)
#fast_design.minimize_bond_lengths(True)

# apply fast design
fast_design.apply(fa_working)

# clone output
design_output = fa_working.clone()

##################################################################
##################################################################
##################################################################
##################################################################
############################# DOCK ###############################
##################################################################
##################################################################
##################################################################
##################################################################

##################################################################
########################### IMPORT ###############################
##################################################################

# load a protein:
designed_pose = fa_working.clone()
native_pose = fa_starting.clone()


##################################################################
################### SWITCH TO CENTROID ###########################
##################################################################

# switch to centroid
switch = SwitchResidueTypeSetMover("centroid")
switch.apply(designed_pose)
switch.apply(native_pose)



##################################################################
################### CHANGE FOLD TREE #############################
##################################################################

# print pdb pose_tree
print(designed_pose.fold_tree())

# change fold tree to keep DEF fixed and move H
setup_foldtree(designed_pose, "DEF_H", Vector1([1])) # "DEF_H"
print(designed_pose.fold_tree())

# set jump information
jump_num = 1

##################################################################
#################### CENTROID ENERGY FUNCTION ####################
##################################################################

# this is the table used for centroid mode docking
motif_dock_score = create_score_function("motif_dock_score")

# what is the score?
print(motif_dock_score(designed_pose)) # zero

##################################################################
############### BROAD LOCAL CENTROID MODE ########################
##################################################################

trans_mag = 3
rot_mag = 8
import rosetta.protocols.rigid as rigid_moves
rigid_body_mover = rigid_moves.RigidBodyPerturbMover(1, rot_mag, trans_mag)
rigid_body_mover.apply(designed_pose)

kT = 0.8
mc = MonteCarlo(designed_pose, motif_dock_score, kT)

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
        trial_mover.apply(designed_pose)
        #pymol.apply(working_pose)

        # store scores and Lrmsd
        scores = np.append(scores, motif_dock_score(designed_pose))
        Lrmsd = np.append(Lrmsd, calc_Lrmsd(designed_pose, native_pose, Vector1([1])))

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
mc.recover_low(designed_pose)

##################################################################
##################### SWITCH BACK TO FULL-ATOM ###################
##################################################################

switch_back = SwitchResidueTypeSetMover("fa_standard")
switch_back.apply(designed_pose)
#pymol.pymol_name("fa_pose")
#pymol.apply(working_pose)


##################################################################
############### FULL-ATOM ENERGY FUNCTION ########################
##################################################################

fa_sfxn = create_score_function("ref2015")


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
mc_fa = MonteCarlo(designed_pose, fa_sfxn, kT)

# trial mover
trial_mover = TrialMover(rigid_body_mover, mc_fa)


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
        trial_mover.apply(designed_pose)
        #pymol.apply(working_pose)

        # store scores and Lrmsd
        scores = np.append(scores, motif_dock_score(designed_pose))
        Lrmsd = np.append(Lrmsd, calc_Lrmsd(designed_pose, native_pose, Vector1([1])))

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
    mc_fa.recover_low(designed_pose)

    # pack the rotamers at the interface
    interface_packer(designed_pose)

mc_fa.show_scores()
mc_fa.recover_low(designed_pose)
#pymol.apply(working_pose)


##################################################################
############### MINIMIZE #########################################
##################################################################
# quick relax...

# Create a MinMover Object to minimize pose
from rosetta.protocols.minimization_packing import *
min_mover = MinMover()
mm_dock = MoveMap()
mm_dock.set_bb(True)
min_mover.movemap(mm_dock)
min_mover.score_function(fa_sfxn)

min_mover.apply(designed_pose)


# find score of best protein
print('best energy: ')
print(fa_sfxn(lowest_energy_pose))
