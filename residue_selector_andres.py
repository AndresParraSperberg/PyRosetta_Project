# import pyrosetta
from pyrosetta import *
init()

###########################################################
#User defined functions go here. Use two spaces after each function.
def print_vector_residues(bool_vec):
    """ Loop through the boolean vector entries to convert True values
        to user-friendly residue numbers. """
    residue_list = []
    for i in range(1, len(bool_vec)+1):
        if bool_vec[i] == 1:
            residue_list.append(i)

    print("The residues belonging to the set are: {}".format(residue_list))


###########################################################
# import and clean pdb (superantigen + TCR + MHC)
from pyrosetta.toolbox import cleanATOM
cleanATOM("2ICW_MHC.pdb")
pose = pose_from_pdb("2ICW_MHC.clean.pdb")

############################################################
#the starting pdb pose
fa_starting = Pose()
fa_starting.assign(pose)

#clone for the analyses
fa_working = fa_starting.clone()

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

############################################################
# setup the packer-task
tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

# These are pretty standard
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())


# Disable repacking/design on non-interface
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
    pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(), not_antigen_helix_neighbors))

# Disable design (just repack) on MHC
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
    pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), combine_DF_selectors))


# Convert the task factory into a PackerTask
packer_task = tf.create_task_and_apply_taskoperations(fa_working)
# View the PackerTask
print(packer_task)

############################################################
# setup the minimizer
mm = pyrosetta.rosetta.core.kinematics.MoveMap()
mm.set_bb(True)
mm.set_chi(True)
mm.set_jump(True)

############################################################
# setup fast design
fast_design = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign(scorefxn_in=scorefxn, standard_repeats=1)
fast_design.cartesian(True)
fast_design.set_task_factory(tf)
fast_design.set_movemap(mm)
fast_design.minimize_bond_angles(True)
fast_design.minimize_bond_lengths(True)

# misc....
loopThrough = pyrosetta.rosetta.core.select.get_residues_from_subset(combine_DF_selectors.apply(fa_working))
for index in loopThrough:
  print(fa_working.pdb_info().pose2pdb(index))




ed
""" We need a section for loop design to patch the missing loop residues
    on the surface of the antigen """

""" We may add a section to do alanine scanning before the protein residues
    are subjected to the design algorithm """


""" We neto decide on what design method to use and what weights/parameters
    should go with that design method. One particular consideration is alanine
    overrepresentation we encountered yesterday. """

""" We will need help regarding how many runs of design/relax/repacking will be
    necessary for a work of this scale. We will also need to find good selection
    criteria to accept/reject the resulting mutants. """

""" We may add a protein-protein docking section for result validation. I'm
    still not convinced that there's no streamlined method for protein docking,
    otherwise why would people even use Rosetta? """

""" We may need a series of functions to generate user-friendly and clear
    output files at the different steps of the calculations. """

""" There may be experimental mutation data available for comparison with our
    results if you feel particularly adventurous. """


'''



# Disable design
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
    pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), not_neighbors))

# Enable design on chain_A_cys_res
# THIS IS WHAT REPLACES THE RESFILE...
aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
aa_to_design.aas_to_keep("ADEFGHIKLMNPQRSTVWY")
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
    aa_to_design, chain_A_cys_res))

'''
