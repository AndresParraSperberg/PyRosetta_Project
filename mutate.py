# to do: minimize full structure, and relax locally with cartesian

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

# relax pose, store lowest energy
# (i.e. import fast relax method from lecture7)
# task factory, what to design, where and how...
# (dock mcm protocol, get interface residues)
# interface energy optimization

##########################
# CreateGlycanSequonMover
# Glycosylate mover
# man5
# glycan tree modeler
# sugarcoat
##########################


from fast_relax_function import *
fast_relax(fa_working, fast_relax_rounds=1)

# find docking score of new mutated protein
print('mutated pose energy: ')
print(scorefxn(fa_working))

