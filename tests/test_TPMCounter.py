#!/usr/bin/env python

import unittest
import logging

from MSTS.Parser.GffGeneParser import GffGeneParser
from MSTS.TPMCounter import TPMCounter

class TestTPMCounter(unittest.TestCase):

    def setUp(self):
        self.inputBam =  "../test-data/input.bam"
        self.inputGff = "../test-data/sample.gff3" 
        self.iTPMCounter = TPMCounter(self.inputBam, self.inputGff)

    def tearDown(self):
        pass

    def test_getGenomeIndex(self):
        """Test getGenomeIndex"""

        lGenes = GffGeneParser(self.inputGff).getAllGenes()
        dExpectedCDSSeq = {}
        dGetSeq = self.iTPMCounter.getGenomeIndex(lGenes,featType="CDS")
        self.assertEquals(dExpectedCDSSeq,dGetSeq)
        dExpectedExonSeq = {}
        dGetSeq = self.iTPMCounter.getGenomeIndex(lGenes,featType="EXON")
        self.assertEquals(dExpectedExonSeq,dGetSeq)

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestTPMCounter)
    unittest.TextTestRunner(verbosity=2).run(suite)
