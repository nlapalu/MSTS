import unittest
import logging

from distutils.cmd import Command
from distutils.core import setup, Extension
from unittest import TextTestRunner, TestLoader, TestSuite

from tests.test_BamConverter import TestBamConverter
from tests.test_ExpNucFileParser import TestExpNucFileParser
from tests.test_FeatureDB import TestFeatureDB
from tests.test_GffGeneParser import TestGffGeneParser
from tests.test_SimpleGffParser import TestSimpleGffParser
from tests.test_TPMCounter import TestTPMCounter


class TestSuite(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):

        logging.disable(logging.CRITICAL)

        lTestCases = [TestBamConverter,TestExpNucFileParser,
                      TestFeatureDB,TestGffGeneParser, TestSimpleGffParser,
                      TestTPMCounter,]

        for case in lTestCases:
            suite = unittest.TestLoader().loadTestsFromTestCase(case)
            t = TextTestRunner(verbosity = 2)
            t.run(suite)
# python2.7
#execfile('MSTS/version.py')
# python 3
exec(open("MSTS/version.py").read())

module = Extension('CMSTS', sources = ["CMSTS.c"])

setup(name='MSTS',
      version = __version__,
      description='MAINE-Seq Tool Suite',
      url='http://github.com/nlapalu/MSTS',
      author='Nicolas Lapalu',
      author_email='nicolas.lapalu@inra.fr',
      license='To define',
      scripts=['bin/MSTS_phasogram.py',
               'bin/MSTS_converter.py',
               'bin/MSTS_feature_phasogram.py',
               'bin/MSTS_count_TPM.py',
               'bin/MSTS_dinuc_frequency.py',
               'bin/MSTS_detect_nucleosomes.py',
               'bin/MSTS_merge_phasograms.py'],
      packages=['MSTS','MSTS.Parser','MSTS.Db','MSTS.Entities',
                'tests','extlib'],
      data_files=[('test-data', ['test-data/input_min.bam',]),
                  ('.',['README.md'])],
      cmdclass={
          'test': TestSuite
      },
      ext_modules = [module]
      )
