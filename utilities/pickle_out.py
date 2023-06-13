import warnings
from m3gnet.models import Relaxer
import ase.io
from ase import Atoms, units
from ase.optimize import BFGS, LBFGS, BFGSLineSearch, QuasiNewton, FIRE
from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG
from ase.constraints import StrainFilter, UnitCellFilter
from ase.io.trajectory import Trajectory
from m3gnet.models._base import Potential
import pickle

element_dat = \
[[ 0  , "X" , 1.0],
[ 1  , "H"  , 0.37],  
[ 2  , "He" , 0.32],  
[ 3  , "Li" , 1.34],  
[ 4  , "Be" , 0.90],  
[ 5  , "B"  , 0.82],  
[ 6  , "C"  , 0.77],  
[ 7  , "N"  , 0.75],  
[ 8  , "O"  , 0.73],  
[ 9  , "F"  , 0.71],  
[ 10 , "Ne" , 0.69],  
[ 11 , "Na" , 1.54],  
[ 12 , "Mg" , 1.30],  
[ 13 , "Al" , 1.18],  
[ 14 , "Si" , 1.11],  
[ 15 , "P"  , 1.06],  
[ 16 , "S"  , 1.02],  
[ 17 , "Cl" , 0.99],  
[ 18 , "Ar" , 0.97],  
[ 19 , "K"  , 1.96],  
[ 20 , "Ca" , 1.74],  
[ 21 , "Sc" , 1.44],  
[ 22 , "Ti" , 1.36],  
[ 23 , "V"  , 1.25],  
[ 24 , "Cr" , 1.27],  
[ 25 , "Mn" , 1.39],  
[ 26 , "Fe" , 1.25],  
[ 27 , "Co" , 1.26],  
[ 28 , "Ni" , 1.21],  
[ 29 , "Cu" , 1.38],  
[ 30 , "Zn" , 1.31],  
[ 31 , "Ga" , 1.26],  
[ 32 , "Ge" , 1.22],  
[ 33 , "As" , 1.19],  
[ 34 , "Se" , 1.16],  
[ 35 , "Br" , 1.14],  
[ 36 , "Kr" , 1.10],  
[ 37 , "Rb" , 2.11],  
[ 38 , "Sr" , 1.92],  
[ 39 , "Y"  , 1.62],  
[ 40 , "Zr" , 1.48],  
[ 41 , "Nb" , 1.37],  
[ 42 , "Mo" , 1.45],  
[ 43 , "Tc" , 1.56],  
[ 44 , "Ru" , 1.26],  
[ 45 , "Rh" , 1.35],  
[ 46 , "Pd" , 1.31],  
[ 47 , "Ag" , 1.53],  
[ 48 , "Cd" , 1.48],  
[ 49 , "In" , 1.44],  
[ 50 , "Sn" , 1.41],  
[ 51 , "Sb" , 1.38],  
[ 52 , "Te" , 1.35],  
[ 53 , "I"  , 1.33],  
[ 54 , "Xe" , 1.30],  
[ 55 , "Cs" , 2.25],  
[ 56 , "Ba" , 1.98],  
[ 57 , "La" , 1.80],  
[ 58 , "Ce" , 1.63],  
[ 59 , "Pr" , 1.76],  
[ 60 , "Nd" , 1.74],  
[ 61 , "Pm" , 1.73],  
[ 62 , "Sm" , 1.72],  
[ 63 , "Eu" , 1.68],  
[ 64 , "Gd" , 1.69],  
[ 56 , "Tb" , 1.68],  
[ 66 , "Dy" , 1.67],  
[ 67 , "Ho" , 1.66],  
[ 68 , "Er" , 1.65],  
[ 69 , "Tm" , 1.64],  
[ 70 , "Yb" , 1.70],  
[ 71 , "Lu" , 1.60],  
[ 72 , "Hf" , 1.50],  
[ 73 , "Ta" , 1.38],  
[ 74 , "W"  , 1.46],  
[ 75 , "Re" , 1.59],  
[ 76 , "Os" , 1.28],  
[ 77 , "Ir" , 1.37],  
[ 78 , "Pt" , 1.28],  
[ 79 , "Au" , 1.44],  
[ 80 , "Hg" , 1.49],  
[ 81 , "Tl" , 1.48],  
[ 82 , "Pb" , 1.47],  
[ 83 , "Bi" , 1.46],  
[ 84 , "Po" , 1.45],  
[ 85 , "At" , 1.47],  
[ 86 , "Rn" , 1.42],  
[ 87 , "Fr" , 2.23],  
[ 88 , "Ra" , 2.01],  
[ 89 , "Ac" , 1.86],  
[ 90 , "Th" , 1.75],  
[ 91 , "Pa" , 1.69],  
[ 92 , "U"  , 1.70],  
[ 93 , "Np" , 1.71],  
[ 94 , "Pu" , 1.72],  
[ 95 , "Am" , 1.66],  
[ 96 , "Cm" , 1.66],  
[ 97 , "Bk" , 1.68],  
[ 98 , "Cf" , 1.68],  
[ 99 , "Es" , 1.65],  
[ 100, "Fm" , 1.67],  
[ 101, "Md" , 1.73],  
[ 102, "No" , 1.76],  
[ 103, "Lr" , 1.61],  
[ 104, "Rf" , 1.57],  
[ 105, "Db" , 1.49],  
[ 106, "Sg" , 1.43],  
[ 107, "Bh" , 1.41],  
[ 108, "Hs" , 1.34],  
[ 109, "Mt" , 1.29],  
[ 110, "Ds" , 1.28],  
[ 111, "Rg" , 1.21],  
[ 112, "Cn" , 1.22]]

results = pickle.load( open( "opt.traj", "rb" ) )
# print("optimized cell=\n", results["cell"][-1])
# print("atomic number=\n", results["atomic_number"])
# print("optimized atom position=\n", results["atom_positions"][-1])
atoms_str = ''
for i in range(len(results["atomic_number"])):
    for j in range(len(element_dat)):
        if int(results["atomic_number"][i]) == int(element_dat[j][0]):
            atoms_str = atoms_str + str(element_dat[j][1])
            
# print("atoms_str=", atoms_str)
atoms = Atoms(atoms_str)
atoms.set_cell(results["cell"][-1])
atoms.set_positions(results["atom_positions"][-1], apply_constraint = False)
atoms.set_atomic_numbers(results["atomic_number"])


ase.io.write('opt.vasp', atoms, direct = True, long_format = True, vasp5 = True)
    