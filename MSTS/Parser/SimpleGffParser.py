#!/usr/bin/env python

import logging
import re
import sys
from sets import Set

from MSTS.Entities.Feature import Feature

class SimpleGffParser(object):

    def __init__(self, inputGffFile="", logLevel='ERROR'):
        """Constructor"""

        self.inputGffFile = inputGffFile
        self.logLevel = logLevel
        logging.basicConfig(level=self.logLevel) 
        
        self.nbFeatures = 0
        self.sReferences = Set()

        try:
            self.filehandle = open(self.inputGffFile, 'r')
        except Exception as e:
            logging.error(e.message)
            sys.exit(1)

    def __del__(self):
        """..."""

        pass

    def _stringToDict(self, string):
        """convert field 9 from string to dict"""

        dAttributes = {}
        for att in string.split(";"):
            tag,val = att.split("=")
            dAttributes[tag] = val.split(",")
        return dAttributes

    def parse(self):
        """Iterator on feature"""

        currentLine = None
        for idx, line in enumerate(self.filehandle):
            currentLine = line.strip()
            if not currentLine:
                pass
            elif re.match('^#', currentLine):
                pass
            else:
		m = re.search(r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+([+-.])\s+(\S+)\s+(\S.*)$", currentLine)
		if m == None:
			raise Exception("Error GFF format line:{}".format(idx))

                id = self._getFeatureTagValue('ID',m.group(9))
	        dAttributes = self._stringToDict(m.group(9))
                score = None
                try:
                    score = float(m.group(6))
                except:
                    pass
                strand = self._getStrand(m.group(7))
                frame = None
                if m.group(8).isdigit():
                    frame = int(m.group(8))
                self.sReferences.add(m.group(1))
                f = Feature(id,m.group(1),m.group(2),m.group(3),int(m.group(4)),int(m.group(5)),score,strand, frame, dAttributes)
                self.nbFeatures += 1

                yield f


    def _getFeatureTagValue(self, tag, line):
        """Return the fist value of the tag property"""
        m = re.search(r".*{mytag}=([^;]*);{{0,1}}.*".format(mytag = tag),line)
        if m:
            return m.group(1).split(',')[0]
        else:
            raise Exception('Cannot find tag {} in string \'{}\''.format(tag, line))


    def _getStrand(self, strand):
        """Return strand as integer(1,-1) instead of +,- """

        if strand == '+':
            return 1
        elif strand == '-':
            return -1
        elif strand == '.':
            return None
        else:
            raise Exception('Cannot defined strand for feature')

    def getsReferences(self):
        """getter set of references"""

        return self.sReferences
