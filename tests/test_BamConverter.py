#!/usr/bin/env python

import unittest
import logging
import sys
import pysam
import array

from MSTS.BamConverter import BamConverter

class TestBamConverter(unittest.TestCase):

    def setUp(self): 

        #self.maxDiff = None
        self.BamFile =  "../test-data/input_min.bam"

    def tearDown(self):
        pass 

    def test_convertWithSingleMode(self):
        """Test convertWithSingleMode"""

        expCurrentSeq = array.array('i',[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 11, 12, 12, 12, 12, 12, 13, 13, 14])
        seq = 'lm_SuperContig_13_v2'
        currentSeq = array.array('i')
        for i in range(0,1634580):
            currentSeq.append(0)
        bc = BamConverter(self.BamFile)
        testCurrentSeq, dFragmentSize, lBedTracks  = bc.convertWithSingleMode(seq, currentSeq)
        self.assertEquals(expCurrentSeq,testCurrentSeq[0:90])


    def test_convertWithFragmentModeFull(self):
        """Test convertWithFragmentModeFull"""

        expCurrentSeq = array.array('i',[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 11, 12, 12, 12, 12, 12, 13, 13, 14])
        seq = 'lm_SuperContig_13_v2'
        currentSeq = array.array('i')
        for i in range(0,1634580):
            currentSeq.append(0)
        bc = BamConverter(self.BamFile, mode="fragment")
        testCurrentSeq, dFragmentSize, lBedTracks  = bc.convertWithFragmentMode(seq, currentSeq)
        self.assertEquals(expCurrentSeq,testCurrentSeq[0:90])
       

    def test_convertWithFragmentModeMiddle(self):
        """Test convertWithFragmentModeMiddle"""

        expCurrentSeq = array.array('i',[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 3, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 1, 1, 0, 0, 3, 0, 2, 2, 1, 0, 1, 0, 1, 0, 1, 3, 1, 0, 1, 0, 1, 1, 0, 1, 2, 0, 0, 4, 2])
        seq = 'lm_SuperContig_13_v2'
        currentSeq = array.array('i')
        for i in range(0,1634580):
            currentSeq.append(0)
        bc = BamConverter(self.BamFile, mode="fragment-middle")
        testCurrentSeq, dFragmentSize, lBedTracks  = bc.convertWithFragmentMode(seq, currentSeq)
        self.assertEquals(expCurrentSeq,testCurrentSeq[0:200])


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBamConverter)
    unittest.TextTestRunner(verbosity=2).run(suite)
