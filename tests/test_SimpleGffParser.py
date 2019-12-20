#!/usr/bin/env python

import unittest

from MSTS.Parser.SimpleGffParser import SimpleGffParser

from MSTS.Entities.Feature import Feature

class TestSimpleGffParser(unittest.TestCase):

    def setUp(self):
        self.parser = SimpleGffParser('test-data/test.simple.gff3', logLevel='DEBUG')
        self.maxDiff = None

    def tearDown(self):
        pass

    def test_parse(self):
        """Test parse"""

        f1 = Feature('gene1','chr1','bioinfobioger','gene',1,200,0.12,1,None,{'toto': ['2', '7'], 'ID': ['gene1'], 'Name': ['gene1']})
        f2 = Feature('gene2','chr1','bioinfobioger','gene',1899,2000,400.0,-1,None,{'ID': ['gene2'], 'Name': ['gene2']})
        f3 = Feature('gene3','chr2','bioinfobioger','gene',18990,20000,40.0,-1,None,{'ID': ['gene3'], 'Name': ['gene3']})
        lExpFeatures = [f1,f2,f3]

        lTestFeatures = []
        for feat in self.parser.parse():
            lTestFeatures.append(feat)
        self.assertEquals(lExpFeatures,lTestFeatures)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSimpleGffParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
