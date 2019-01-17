# import pyrosetta
from pyrosetta import *
init()

def print_vector_residues(bool_vec):
    residue_list = []
    for i in range(1, len(bool_vec)+1):
        if bool_vec[i] == 1:
            residue_list.append(i)
            
    print("The residues belonging to the set are: {}".format(residue_list))



# import and clean pdb (superantigen + TCR + MHC)
from pyrosetta.toolbox import cleanATOM
cleanATOM("2ICW_MHC.pdb")
pose = pose_from_pdb("2ICW_MHC.clean.pdb")

############################################################
#the starting pdb pose
fa_starting = Pose()
fa_starting.assign(pose)

#cone for the analyses
fa_working = fa_starting.clone()

#find the pose numbering from pdb numbering
start_residue = fa_working.pdb_info().pdb2pose('D', 55)
end_residue = fa_working.pdb_info().pdb2pose('D', 72)

#first define the main helix of the antigen in contact with the protein
antigen_helix_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
antigen_helix_selector.set_index("{}-{}".format(start_residue, end_residue))
antigen_helix_residue_vector = antigen_helix_selector.apply(fa_working)
#print_vector_residues(antigen_helix_residue_vector)

#find the neighbors
antigen_helix_neighbor_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(antigen_helix_residue_vector, 8, True)
antigen_helix_neighbors = antigen_helix_neighbor_selector.apply(fa_working)

#define the chain H residues belonging to the protein
H_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('H')

combine_selectors = pyrosetta.rosetta.core.select.residue_selector.AND_combine(antigen_helix_neighbor_selector, H_selector)
H_residues = combine_selectors.apply(fa_working)
print_vector_residues(H_residues)
#will use and_combine to chose within range and in chain h
