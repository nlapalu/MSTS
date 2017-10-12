import unittest
import logging

from distutils.cmd import Command
from distutils.core import setup, Extension
from unittest import TextTestRunner, TestLoader, TestSuite

#from tests.test_AlignDB import TestAlignDB


#try:
#    from tests.test_BlastXMLParser import TestBlastXMLParser
#    Biopython_available = True
#except ImportError:
#    Biopython_available = False

class TestSuite(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):

        logging.disable(logging.CRITICAL) 

        lTestCases = [TestAlignDB,]

#        if Biopython_available:
#            lTestCases.append(TestBlastXMLParser)
        for case in lTestCases:
            suite = unittest.TestLoader().loadTestsFromTestCase(case)
            t = TextTestRunner(verbosity = 2)
            t.run(suite)

execfile('MSTS/version.py')

module = Extension('CMSTS', sources = ["CMSTS.c"])

setup(name='MSTS',
      version = __version__,
      description='MAINE-Seq Tool Suite',
      url='http://github.com/nlapalu/MSTS',
      author='Nicolas Lapalu',
      author_email='nicolas.lapalu@inra.fr',
      license='To define',
      scripts=['bin/MSTS_phasogram.py',
               'bin/MSTS_converter.py'],
      packages=['MSTS','tests',],
      data_files=[('test-data', ['test-data/input_min.bam',]),
                  ('.',['README.md'])],
      cmdclass={
          'test': TestSuite
      },
      ext_modules = [module]
      )
