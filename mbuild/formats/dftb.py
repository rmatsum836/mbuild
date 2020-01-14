import numpy as np
import json
import mbuild as mb

def read_dftb(hsd_file, compound=None):
    """
    Read the atom information from a DFTB+ HSD File

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

    return compound
