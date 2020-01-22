import numpy as np
import json
import mbuild as mb
from decimal import Decimal

def read_dftb(gen_file, compound=None):
    """
    Read the atom information from a DFTB+ Gen File
    Currently only reading in 'GenFormat' types

    Parameters
    ----------
    gen_file: str
        Input file for DFTB+
    """

    if compound is None:
        compound = mb.Compound()

    breaks = False
    with open(gen_file, 'r') as fi:
        line = fi.readline()
        # Read in number of atoms and type of geometry
        n_atoms = int(line.split()[0])
        geo_type = line.split()[1]
        # Read in only the geometry section of the file
        geo_info = fi.readline().split()
        # Read in geometry coordinates
        geo_coords = fi.readlines()[:n_atoms]
        # Map atom to a type (integer)
        atomtype_dict = {i+1:j for i,j in zip(range(len(geo_info)), geo_info)}
    
    coords = np.zeros(shape=(n_atoms, 3), dtype=np.float64)
    if geo_type == 'C':
        for row, line in enumerate(geo_coords):
            # Grab coordinates and convert to floats
            # Convert from angstroms to nm
            coords[row] = [float(i)/10 for i in line.split()[-3:]]

            # Add particle
            particle = mb.Compound(pos=coords[row], name=atomtype_dict[int(line.split()[1])])
            compound.add(particle)

    # TODO: Read in supercell lattice information

    return compound

def write_dftb(structure, filename, geometry='C'):
    """
    Write out ParmEd structure to a DFTB+ GEN File

    Parameters
    ----------
    structure: ParmEd Structure
    filename: str
        Output file for DFTB+

    geometry: str, default='C'
        Type of geometry:
        'C' is cluster (non-periodic)
        'F' is supercell in fractions of lattice vectors
        'S' is for supercell in Cartesian Coordinates
    """
    if isinstance(structure, mb.Compound):
        structure = structure.to_parmed()

    xyz = np.array([[atom.xx*10, atom.xy*10, atom.xz*10]
        for atom in structure.atoms])

    # TODO: Convert cartesian coordinates
    if geometry in ['F', 'S']:
        pass
    names = [i.name for i in structure.atoms]
    n_atoms = len(structure.atoms)
    elements = list(set([i.name for i in structure.atoms]))

    element_dict = {i:j+1 for i,j in zip(elements, range(len(elements)))}

    with open(filename, 'w') as data:
        data.write('{} {}\n'.format(n_atoms, geometry))
        for element in element_dict.keys():
            data.write('{} '.format(element))
        data.write('\n')
        for idx in range(len(names)):
            data.write('{} {}  {:.11E} {:.11E} {:.11E}\n'.format(
                idx+1, element_dict[names[idx]], Decimal(xyz[idx][0]),
                Decimal(xyz[idx][1]), Decimal(xyz[idx][2])))

    # TODO: Write supercell info
