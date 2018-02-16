#!/usr/bin/env python

import unittest
import logging

from MSTS.Entities.Gene import Gene

class TestGene(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_getExonTotalLength(self):
        """Test getExonTotalLength"""


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestGene)
    unittest.TextTestRunner(verbosity=2).run(suite)
