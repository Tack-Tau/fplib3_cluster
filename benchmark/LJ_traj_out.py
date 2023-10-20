import os
import sys
import numpy as np
import ase.io
from ase.io.trajectory import Trajectory
from ase.calculators.lj import LennardJones

fp_traj = Trajectory('fp_opt.traj')
LJ_traj = Trajectory('opt.traj')
fp_count = 0
LJ_count = 0
fp_max = 0
LJ_max = 0
for atoms in fp_traj:
    fp_max = fp_max + 1

for atoms in LJ_traj:
    LJ_max = LJ_max + 1

max_count = max(fp_max, LJ_max)

for atoms in fp_traj:
    fp_count = fp_count + 1
    if fp_count <= max_count:
        f_max = np.amax( np.absolute( atoms.get_forces() ) )
        calc = LennardJones()
        calc.parameters.epsilon = 1.0
        calc.parameters.sigma = 1.0
        calc.parameters.rc = 1000.0
        calc.parameters.smooth = False
        atoms.calc = calc
        e = atoms.get_potential_energy() / len(atoms)
        # f_max = np.amax( np.absolute( atoms.get_forces() ) )
        print(fp_count, e, f_max)
    
for atoms in LJ_traj:
    LJ_count = LJ_count + 1
    if LJ_count <= max_count:
        e = atoms.get_potential_energy() / len(atoms)
        f_max = np.amax( np.absolute( atoms.get_forces() ) )
        print(LJ_count, e, f_max)