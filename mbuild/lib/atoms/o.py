import numpy as np

import mbuild as mb


class O(mb.Compound):
    """A hydrogen atom with two overlayed ports."""
    def __init__(self):
        super(O, self).__init__()
        self.add(mb.Particle(name='O'))

        self.add(mb.Port(anchor=self[0]), 'up')
        self.add(mb.Port(anchor=self[0]), 'down')
        self['up'].spin(np.pi*3/4, [0, 0, 1])
        self['up'].translate(np.array([0, 0.07, 0]))
        self['down'].spin(-np.pi*3/4, [0, 0, 1])
        self['down'].translate([0, 0.07, 0])


if __name__ == '__main__':
    m = H()
    m.save('o.mol2', overwrite=True)
