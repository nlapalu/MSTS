#!/usr/bin/env python

import unittest
import logging
import sys
import os

from MSTS.Parser.ExpNucFileParser import ExpNucFileParser

class TestExpNucFileParser(unittest.TestCase):

    def setUp(self):
        self.parser = ExpNucFileParser('test-data/expnucfile.txt', logLevel='DEBUG')
        self.maxDiff = None

    def tearDown(self):
        pass

    def test_parse(self):
        """test parse function"""

        dParsedConditions = self.parser.parse()
        dExpectedConditons = {'c1': [('1', 'nucfile1.1', 'bwfile1.1'), 
                                     ('2', 'nucfile1.2', 'bwfile1.2')],
                              'c2': [('3', 'nucfile2.1', 'bwfile2.1'),
                                     ('4', 'nucfile2.2', 'bwfile2.2'),
                                     ('5', 'nucfile2.3', 'bwfile2.3')]}
        self.assertEquals(dExpectedConditons,dParsedConditions)

    def test_getlLabels(self):
        """test getlLabels"""

        lExpectedLabels = ['1','2','3','4','5']
        dParsedConditions = self.parser.parse()
        lParsedLabels = self.parser.getlLabels()
        self.assertEquals(lExpectedLabels,lParsedLabels)

    def test_getlConditions(self):
        """test getlConditions"""

        lExpectedlConditions = [[0,1],[2,3,4]]

        self.parser.parse()
        lParsedlConditions = self.parser.getlConditions()
        self.assertEquals(lExpectedlConditions,lParsedlConditions)
        
    def test_getlNucFiles(self):
        """test getlNucFiles"""

        lExpectedlNucFiles = ['nucfile1.1','nucfile1.2','nucfile2.1','nucfile2.2','nucfile2.3']
        self.parser.parse()
        lParsedlNucFiles = self.parser.getlNucFiles()
        self.assertEquals(lExpectedlNucFiles,lParsedlNucFiles)

    def test_getlWigFiles(self):
        """test getlWigFiles"""

        lExpectedlWigFiles = ['bwfile1.1','bwfile1.2','bwfile2.1','bwfile2.2','bwfile2.3']
        self.parser.parse()
        lParsedlWigFiles = self.parser.getlWigFiles()
        self.assertEquals(lExpectedlWigFiles,lParsedlWigFiles)




if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestExpNucFileParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
