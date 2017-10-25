#!/usr/bin/env python

import unittest

from MSTS.Db.FeatureDB import FeatureDB

from MSTS.Entities.Feature import Feature

class TestFeatureDB(unittest.TestCase):

    def setUp(self):
        self.db = FeatureDB(dbfile='test.db', logLevel='DEBUG')
        self.maxDiff = None

    def tearDown(self):
#        pass
        self.db.deleteDB()

    def test_selectAllFeatures(self):
        """Test selectAllGenes"""

        lFeatures = []
        feature = Feature('G00001','Chr1','bioinfobioger','gene',23988,24919,50.0,-1,None,{'ID':'G00001','Name':'G00001'})
        lFeatures.append(feature)
        self.db.insertlFeatures(lFeatures)

        self.assertEquals(lFeatures,self.db.selectAllFeatures())

# TODO : test selectFeatureTypeFromReference


        # self.db.deleteAllGenes()  TODO


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFeatureDB)
    unittest.TextTestRunner(verbosity=2).run(suite)
