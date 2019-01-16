from pyrosetta import *
init()

def fast_relax(pose, fast_relax_rounds=1, use_ex_rotamers=True):
    '''Relax a pose.'''
    rosetta.basic.options.set_boolean_option('relax:constrain_relax_to_start_coords', True)
    rosetta.basic.options.set_boolean_option('relax:ramp_constraints', True)

    task_factory = rosetta.core.pack.task.TaskFactory()

    if use_ex_rotamers:
        ers = rosetta.core.pack.task.operation.ExtraRotamersGeneric()
        ers.ex1(True)
        ers.ex2(True)
        ers.extrachi_cutoff(18)
        task_factory.push_back(ers)

    task_factory.push_back(rosetta.core.pack.task.operation.RestrictToRepacking())
    lac = rosetta.protocols.task_operations.LimitAromaChi2Operation()
    task_factory.push_back(lac)

    sfxn = rosetta.core.scoring.get_score_function()
    fast_relax_mover = rosetta.protocols.relax.FastRelax(sfxn, fast_relax_rounds)

    fast_relax_mover.apply(pose)


