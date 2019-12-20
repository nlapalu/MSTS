#!/usr/bin/env python

import unittest
import logging
import sys
import os

from MSTS.Parser.GffGeneParser import GffGeneParser
from MSTS.Entities.Gene import Gene
from MSTS.Entities.Transcript import Transcript
from MSTS.Entities.CDS import CDS
from MSTS.Entities.Exon import Exon

class TestGffGeneParser(unittest.TestCase):

    def setUp(self):
        self.parser = GffGeneParser('test-data/sample.gff3', logLevel='DEBUG')
        self.maxDiff = None

    def tearDown(self):
        pass

    def test_getAllGenes(self):
        """test  """

        lExpectedGenes = [Gene("G1","seq1",19966,20759,-1,[Transcript("T1","seq1",19966,20759,-1,"G1",[CDS("G1C1","seq1",20186,20668,-1,"T1")],[Exon("G1E1","seq1",19966,20130,-1,["T1"]),Exon("G1E2","seq1",20186,20759,-1,["T1"])])]),
                          Gene("G4","seq2",139543,140931,1,[Transcript("T4","seq2",139543,140931,1,"G4",[CDS("G4C1","seq2",139543,139555,1,"T4"),CDS("G4C1","seq2",139642,139709,1,"T4"),CDS("G4C1","seq2",139790,140054,1,"T4"),CDS("G4C1","seq2",140150,140931,1,"T4")],[Exon("G4E1","seq2",139543,139555,1,["T4","T4.2"]),Exon("G4E2","seq2",139642,139709,1,["T4"]),Exon("G4E3","seq2",139790,140054,1,["T4"]),Exon("G4E4","seq2",140150,140931,1,["T4"])]),Transcript("T4.2","seq2",139543,143300,1,"G4",[CDS("G4C2","seq2",139543,139555,1,"T4.2"),CDS("G4C2","seq2",139642,139700,1,"T4.2"), CDS("G4C2","seq2",140900,143250,1,"T4.2")],[Exon("G4E1","seq2",139543,139555,1,["T4","T4.2"]),Exon("G4E2.2","seq2",139642,139700,1,["T4.2"]),Exon("G4E4.2","seq2",140900,143300,1,["T4.2"])])])
                         ]
        lParsedGenes = self.parser.getAllGenes()
        self.assertEquals(lExpectedGenes[0],lParsedGenes[0])
        self.assertEquals(lExpectedGenes[1],lParsedGenes[3])

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestGffGeneParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
