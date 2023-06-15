import os
import sys
import numpy as np
import ase.io
from ase.optimize import BFGS, LBFGS, BFGSLineSearch, QuasiNewton, FIRE
from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG
from ase.constraints import StrainFilter, UnitCellFilter
from ase.io.trajectory import Trajectory


atoms = ase.io.read('.'+'/'+'POSCAR')
atoms.set_pbc((False, False, False)) # For clusters turn off PBC
# ase.io.write('input.xyz', atoms, plain = True)
ase.io.vasp.write_vasp('input.vasp', atoms, direct = True)
trajfile = 'opt.traj'
print("Number of atoms:", len(atoms))


from ase.calculators.lj import LennardJones
calc1 = LennardJones()
calc1.parameters.epsilon = 1.0
calc1.parameters.sigma = 1.0
calc1.parameters.rc = 3.0
calc1.parameters.smooth = True

atoms.calc = calc1
print ("LJ_energy:\n", atoms.get_potential_energy())
print ("LJ_forces:\n", atoms.get_forces())
# print ("LJ_stress:\n", atoms.get_stress())

'''
from SF_LJ_api4ase import ShiftedForceLennardJones

calc1 = ShiftedForceLennardJones()
calc1.parameters.epsilon = np.array([1.00, 1.00, 1.00])
calc1.parameters.sigma = np.array([1.00, 1.00, 1.00])
calc1.parameters.rc = 3.0 * np.array([1.00, 1.00, 1.00])

atoms.calc = calc1
print ("SFLJ_energy:\n", atoms.get_potential_energy())
print ("SFLJ_forces:\n", atoms.get_forces())
# print ("SFLJ_stress:\n", atoms.get_stress())
'''

##################################################################################################
# Sigma gives a measurement of how close two nonbonding particles can get and is thus referred to as the van der Waals radius. It is equal to one-half of the internuclear distance between nonbonding particles.
# Ideally, r_min == 2**(1/6) * sigma == 2.0 * r_cov, which means van der Waals radius is approximately two times larger than covalent radius.

# Reference:
# https://en.wikipedia.org/wiki/Lennard-Jones_potential
# https://en.wikipedia.org/wiki/Van_der_Waals_radius
# https://en.wikipedia.org/wiki/Covalent_radius
##################################################################################################


from fplib3_cluster_api4ase import fp_GD_Calculator
from fplib3_cluster_mixing import MixedCalculator

calc2 = fp_GD_Calculator(
            cutoff = 10.0,
            contract = False,
            znucl = np.array([1], int),
            lmax = 0,
            nx = 50, # For clusters choose nx<=len(atoms)
            ntyp = 1
            )

atoms.calc = calc2
print ("fp_energy:\n", atoms.get_potential_energy())
print ("fp_forces:\n", atoms.get_forces())
# print ("fp_stress:\n", atoms.get_stress())


calc = MixedCalculator(calc1, calc2)
atoms.calc = calc
print ("mixed_energy:\n", atoms.get_potential_energy())
print ("mixed_forces:\n", atoms.get_forces())
# print ("mixed_stress:\n", atoms.get_stress())


############################## Relaxation type ##############################
#     https ://wiki.fysik.dtu.dk/ase/ase/optimize.html#module-optimize      #
#     https ://wiki.fysik.dtu.dk/ase/ase/constraints.html                   #
#############################################################################

af = atoms
# af = StrainFilter(atoms)
# af = UnitCellFilter(atoms, scalar_pressure = 0.062415)
# af = UnitCellFilter(atoms, scalar_pressure = 0.0)

############################## Relaxation method ##############################

# opt = BFGS(af, maxstep = 1.e-1, trajectory = trajfile)
opt = FIRE(af, maxstep = 1.e-1, trajectory = trajfile)
# opt = LBFGS(af, maxstep = 1.e-1, trajectory = trajfile, memory = 10, use_line_search = True)
# opt = LBFGS(af, maxstep = 1.e-1, trajectory = trajfile, memory = 10, use_line_search = False)
# opt = SciPyFminCG(af, maxstep = 1.e-1, trajectory = trajfile)
# opt = SciPyFminBFGS(af, maxstep = 1.e-1, trajectory = trajfile)

opt.run(fmax = 1.e-3, steps = 5000)

traj = Trajectory(trajfile)
atoms_final = traj[-1]
# ase.io.write('opt.xyz', atoms_final, plain = True)
ase.io.write('opt.vasp', atoms_final, direct = True, long_format = True, vasp5 = True)

final_structure = atoms.get_positions()
final_energy_per_atom = float( atoms.get_potential_energy() / len(atoms_final) )

print("Relaxed structure in Cartesian coordinates is \n{0:s}".\
      format(np.array_str(final_structure, precision=6, suppress_small=False)))
print("Final energy per atom is \n{0:.6f}".format(final_energy_per_atom))
