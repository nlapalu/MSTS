#!/usr/bin/env python

import unittest
import logging

from MSTS.Parser.GffGeneParser import GffGeneParser
from MSTS.TPMCounter import TPMCounter

class TestTPMCounter(unittest.TestCase):

    def setUp(self):
        self.inputBam =  "test-data/input.bam"
        self.inputGff = "test-data/sample.gff3"
        self.iTPMCounter = TPMCounter(self.inputBam, self.inputGff)

    def tearDown(self):
        pass

    def no_test_getGenomeIndex(self):
        """Test getGenomeIndex"""

        lGenes = GffGeneParser(self.inputGff).getAllGenes()
#        print lGenes
        dExpectedCDSSeq = {}
        dGetSeq = self.iTPMCounter.getGenomeIndex(lGenes,featType="CDS")
        self.assertEquals(dExpectedCDSSeq,dGetSeq)
        dExpectedExonSeq = {}
        dGetSeq = self.iTPMCounter.getGenomeIndex(lGenes,featType="EXON")
        self.assertEquals(dExpectedExonSeq,dGetSeq)

    def test_getTranscriptIndex(self):
        """Test getGenomeIndex"""

        lGenes = GffGeneParser(self.inputGff).getAllGenes()
        for gene in lGenes:
            for transcript in gene.lTranscripts:
                sumExpected = transcript.getLength() - transcript.getExonTotalLength() 
                #print transcript.id, sumExpected
                dGetSeq = self.iTPMCounter.getTranscriptIndex(transcript,featType="EXON")
                sumGet = sum(dGetSeq[transcript.seqid][transcript.start-1:transcript.end])
                self.assertEquals(sumExpected,sumGet)
                sumExpected = transcript.getLength() - transcript.getCDSTotalLength() 
                #print transcript.id, sumExpected
                dGetSeq = self.iTPMCounter.getTranscriptIndex(transcript,featType="CDS")
                sumGet = sum(dGetSeq[transcript.seqid][transcript.start-1:transcript.end])
                self.assertEquals(sumExpected,sumGet)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestTPMCounter)
    unittest.TextTestRunner(verbosity=2).run(suite)
