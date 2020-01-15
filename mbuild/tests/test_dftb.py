import numpy as np
import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn
from mbuild.formats.dftb import read_dftb, write_dftb
from mbuild.utils.io import has_foyer

@pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
class Testdftb(BaseTest):
    def test_save(self, ethane):
        ethane.save(filename='ethane.gen')
