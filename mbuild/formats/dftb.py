import numpy as np
import json
import mbuild as mb
from decimal import Decimal

def read_dftb(hsd_file, compound=None):
    """
    Read the atom information from a DFTB+ HSD File
    Currently only reading in 'GenFormat' types

    Parameters
    ----------
    hsd_file: str
        Input file for DFTB+
    """

    if compound is None:
        compound = mb.Compound()

    breaks = False
    with open(hsd_file, 'r') as fi:
        for i, line in enumerate(fi):
            if breaks == False:
                if line.split()[0] == 'Geometry':
                    breaks = True
            else:
                # TODO: Handle geometry from a file
                if 'gen' in line.split()[0]:
                    continue
                else:
                    n_atoms = int(line.split()[0])
                    break

    # Read in only the geometry section of the file
    # Not sure if skipping 3 lines from the geometry line will work for every file
    geo_info = open(hsd_file, 'r').readlines()[i+1].split()
    geo_coords = open(hsd_file, 'r').readlines()[i+3:i+n_atoms+3]

    # Map atom to a type (integer)
    atomtype_dict = {i+1:j for i,j in zip(range(len(geo_info)), geo_info)}

    coords = np.zeros(shape=(n_atoms, 3), dtype=np.float64)
    for row, line in enumerate(geo_coords):
        # Grab coordinates and convert to floats
        # Convert from angstroms to nm
        coords[row] = [float(i)/10 for i in line.split()[-3:]]

        # Add particle
        particle = mb.Compound(pos=coords[row], name=atomtype_dict[int(line.split()[1])])
        compound.add(particle)

    # TODO: Read in supercell lattice information

    return compound

def write_dftb(structure, filename, geometry='S'):
    """
    Write out ParmEd structure to a DFTB+ GEN File

    Parameters
    ----------
    structure: ParmEd Structure
    filename: str
        Output file for DFTB+

    geometry: str, default='S'
        Type of geometry, 'F' is supercell in fractions of lattice vectors,
        'S' is for supercell in Cartesian Coordinates
    """

    xyz = np.array([[atom.xx*10, atom.xy*10, atom.xz*10]
        for atom in structure.atoms])

    # TODO: Convert cartesian coordinates
    if geometry == 'F':
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
