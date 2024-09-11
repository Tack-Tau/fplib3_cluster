import sys
import os
from ase.io import read, write
from pymatgen.core import Molecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer

def convert_vasp_to_xyz(vasp_file, xyz_file):
    """
    Convert a VASP POSCAR structure to XYZ format with PBC set to False.

    Args:
    vasp_file (str): Path to the VASP POSCAR file (*.vasp).
    xyz_file (str): Path to save the converted XYZ file.
    """
    # Read the VASP file using ASE
    structure = read(vasp_file)

    # Set periodic boundary conditions (PBC) to False
    structure.set_pbc(False)

    # Write the structure to an XYZ file
    write(xyz_file, structure)
    print(f"Converted {vasp_file} to {xyz_file} with PBC set to False.")

def analyze_point_group(xyz_file):
    """
    Analyze the point group symmetry of a molecule using Pymatgen and return the point group and number of symmetry operations.

    Args:
    xyz_file (str): Path to the XYZ file.

    Returns:
    str: The point group symmetry of the molecule.
    int: The number of symmetry operations.
    """
    # Read the XYZ file as a Molecule using Pymatgen
    molecule = Molecule.from_file(xyz_file)

    # Analyze the point group symmetry
    pga = PointGroupAnalyzer(molecule)
    point_group = pga.get_pointgroup()

    # Get the number of symmetry operations
    num_symmetry_operations = len(pga.get_symmetry_operations())

    return point_group, num_symmetry_operations

if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 2:
        print("Usage: python script.py <structure_file.vasp>")
        sys.exit(1)

    # Get the structure file name from the command-line arguments
    vasp_file = sys.argv[1]

    # Check if the file exists
    if not os.path.isfile(vasp_file):
        print(f"Error: The file {vasp_file} does not exist.")
        sys.exit(1)

    # Define the output XYZ file name
    xyz_file = os.path.splitext(vasp_file)[0] + ".xyz"

    # Step 1: Convert VASP to XYZ with PBC set to False
    convert_vasp_to_xyz(vasp_file, xyz_file)

    # Step 2: Analyze the point group symmetry using Pymatgen and output the number of symmetry operations
    point_group, num_symmetry_operations = analyze_point_group(xyz_file)

    # Output the point group symmetry and the number of symmetry operations
    print(f"Point group: {point_group}, {num_symmetry_operations} symmetry operations")

