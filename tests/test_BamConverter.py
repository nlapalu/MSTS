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
        self.BamFile2 =  "../test-data/input_micro.bam"

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

        expCurrentSeq = array.array('i',[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 2, 0, 1, 0, 0, 0, 1, 0, 1, 2, 0, 0, 2, 0, 3, 0, 1, 0, 4, 0, 0, 3, 3, 0, 0, 1, 3, 1, 2, 0, 3, 0, 3, 4, 1, 0, 1, 1, 1, 1, 4, 3, 1, 1, 1, 0, 3, 1, 0, 2, 3, 0, 3, 5])
        seq = 'lm_SuperContig_13_v2'
        currentSeq = array.array('i')
        for i in range(0,1634580):
            currentSeq.append(0)
        bc = BamConverter(self.BamFile, mode="fragment-middle")
        testCurrentSeq, dFragmentSize, lBedTracks  = bc.convertWithFragmentMode(seq, currentSeq)
        self.assertEquals(expCurrentSeq,testCurrentSeq[0:200])

    def test_convertWithFragmentModeMiddleWindow(self):
        """Test convertWithFragmentModeMiddle with window"""

        expCurrentSeq = array.array('i',[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 2, 2, 2, 1, 2, 1, 2, 1, 1, 2, 2, 3, 1, 1, 0, 1, 1, 2, 3, 3, 2, 2, 2, 5, 3, 4, 1, 5, 4, 4, 3, 6, 6, 3, 1, 4, 5, 6, 3, 5, 3, 6, 7, 8, 5, 2, 2, 3, 3, 6, 8, 8, 5, 3, 2, 4, 4, 4, 3, 5, 5, 6, 8, 10])
        seq = 'lm_SuperContig_13_v2'
        currentSeq = array.array('i')
        for i in range(0,1634580):
            currentSeq.append(0)
        bc = BamConverter(self.BamFile, mode="fragment-middle", window=1)
        testCurrentSeq, dFragmentSize, lBedTracks  = bc.convertWithFragmentMode(seq, currentSeq)
        self.assertEquals(expCurrentSeq,testCurrentSeq[0:200])



if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBamConverter)
    unittest.TextTestRunner(verbosity=2).run(suite)
