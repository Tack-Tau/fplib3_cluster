# Fingerprint library for clusters
### This repo serves as an independent nonperiodic system implementation of the original Fingerprint library ([fplib3](https://github.com/Tack-Tau/fplib3))

## Dependencies
* Python >= 3.8.5
* Numpy >= 1.21.4
* Scipy >= 1.8.0
* Numba >= 0.56.2
* ASE >= 3.22.1

## Setup
For conda installation guide please visit their [website](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) \
`conda create -n fplibenv python=3.8 pip ; conda activate fplibenv`\
`python3 -m pip install --user -U pip setuptools wheel numpy scipy ase numba`\
`git clone https://github.com/Tack-Tau/fplib3.git ./fplib3`

## Usage
### Basic ASE style documentation
See details for [ASE calculator class](https://wiki.fysik.dtu.dk/ase/development/calculators.html)
and [ASE calculator proposal](https://wiki.fysik.dtu.dk/ase/development/proposals/calculators.html#aep1)
```
    Fingerprint Calculator interface for ASE
    
        Implemented Properties:
        
            'energy': Sum of atomic fingerprint distance (L2 norm of two atomic 
                                                          fingerprint vectors)
            
            'forces': Gradient of fingerprint energy, using Hellmann–Feynman theorem
            
        Parameters:
        
            atoms:  object
                Attach an atoms object to the calculator.
                
            contract: bool
                Calculate fingerprint vector in contracted Guassian-type orbitals or not
            
            ntype: int
                Number of different types of atoms in unit cell
            
            nx: int
                Maximum number of atoms in the sphere with cutoff radius for specific cell site
                
            lmax: int
                Integer to control whether using s orbitals only or both s and p orbitals for 
                calculating the Guassian overlap matrix (0 for s orbitals only, other integers
                will indicate that using both s and p orbitals)
                
            cutoff: float
                Cutoff radius for f_c(r) (smooth cutoff function) [amp], unit in Angstroms
                For clusters please choose as large as possible to include all the atoms in the structure
                
```


### Calling fplib3 calculator from ASE API
```python
import numpy as np
import ase.io
from ase.optimize import BFGS, LBFGS, BFGSLineSearch, QuasiNewton, FIRE
from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG
from ase.constraints import StrainFilter, UnitCellFilter
from ase.io.trajectory import Trajectory

from fplib3_api4ase import fp_GD_Calculator
# from fplib3_mixing import MixedCalculator
# from ase.calculators.mixing import MixedCalculator
# from ase.calculators.vasp import Vasp

atoms = ase.io.read('.'+'/'+'POSCAR')
ase.io.vasp.write_vasp('input.vasp', atoms, direct=True)
trajfile = 'opt.traj'

from functools import reduce

chem_nums = list(atoms.numbers)
znucl_list = reduce(lambda re, x: re+[x] if x not in re else re, chem_nums, [])
ntyp = len(znucl_list)
znucl = np.array(znucl_list, int)

calc = fp_GD_Calculator(
            cutoff = 10.0,
            contract = False,
            znucl = znucl,
            lmax = 0,
            nx = 50, # For clusters choose nx<=len(atoms)
            ntyp = ntyp
            )

atoms.calc = calc

# calc.test_energy_consistency(atoms = atoms)
# calc.test_force_consistency(atoms = atoms)

print ("fp_energy:\n", atoms.get_potential_energy())
print ("fp_forces:\n", atoms.get_forces())
# print ("fp_stress:\n", atoms.get_stress())          # No stress method implemented

af = atoms
# af = StrainFilter(atoms)                            # Do not use strain filter for clusters
# af = UnitCellFilter(atoms, scalar_pressure = 0.0)   # Do not use unit cell filter for clusters

############################## Relaxation method ##############################\

# opt = BFGS(af, maxstep = 1.e-1, trajectory = trajfile)
opt = FIRE(af, maxstep = 1.e-1, trajectory = trajfile)
# opt = LBFGS(af, maxstep = 1.e-1, trajectory = trajfile, memory = 10, use_line_search = True)
# opt = LBFGS(af, maxstep = 1.e-1, trajectory = trajfile, memory = 10, use_line_search = False)
# opt = SciPyFminCG(af, maxstep = 1.e-1, trajectory = trajfile)
# opt = SciPyFminBFGS(af, maxstep = 1.e-1, trajectory = trajfile)

opt.run(fmax = 1.e-5)

traj = Trajectory(trajfile)
ase.io.write('opt.vasp', traj[-1], direct = True, long_format=True, vasp5 = True)
```
## Citation
If you use this Fingerprint Library (or modified version) for your research please kindly cite our paper:
```
@article{taoAcceleratingStructuralOptimization2024,
  title = {Accelerating Structural Optimization through Fingerprinting Space Integration on the Potential Energy Surface},
  author = {Tao, Shuo and Shao, Xuecheng and Zhu, Li},
  year = {2024},
  month = mar,
  journal = {J. Phys. Chem. Lett.},
  volume = {15},
  number = {11},
  pages = {3185--3190},
  doi = {10.1021/acs.jpclett.4c00275},
  url = {https://pubs.acs.org/doi/10.1021/acs.jpclett.4c00275}
}
```
If you use Fingerprint distance as a metric to measure crystal similarity please also cite the following paper:
```
@article{zhuFingerprintBasedMetric2016,
  title = {A Fingerprint Based Metric for Measuring Similarities of Crystalline Structures},
  author = {Zhu, Li and Amsler, Maximilian and Fuhrer, Tobias and Schaefer, Bastian and Faraji, Somayeh and Rostami, Samare and Ghasemi, S. Alireza and Sadeghi, Ali and Grauzinyte, Migle and Wolverton, Chris and Goedecker, Stefan},
  year = {2016},
  month = jan,
  journal = {The Journal of Chemical Physics},
  volume = {144},
  number = {3},
  pages = {034203},
  doi = {10.1063/1.4940026},
  url = {https://doi.org/10.1063/1.4940026}
}
```
If you use GOM Fingerprint as a descriptor for developing MLP and other machine-learning related research please also cite the following paper:
```
@article{sadeghiMetricsMeasuringDistances2013,
  title = {Metrics for Measuring Distances in Configuration Spaces},
  author = {Sadeghi, Ali and Ghasemi, S. Alireza and Schaefer, Bastian and Mohr, Stephan and Lill, Markus A. and Goedecker, Stefan},
  year = {2013},
  month = nov,
  journal = {The Journal of Chemical Physics},
  volume = {139},
  number = {18},
  pages = {184118},
  doi = {10.1063/1.4828704},
  url = {https://pubs.aip.org/aip/jcp/article/317391}
}
```
