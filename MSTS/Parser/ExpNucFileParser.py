#!/usr/bin/env python

import logging
import re
import sys

class ExpNucFileParser(object):

    def __init__(self, inputFile="", logLevel='ERROR'):
        """Constructor"""

        self.inputFile = inputFile
        self.logLevel = logLevel
        logging.basicConfig(level=self.logLevel) 
        
        self.dConditions = {}

        try:
            self.filehandle = open(self.inputFile, 'r')
        except Exception as e:
            logging.error(e)
            sys.exit(1)

    def parse(self):
        """parse file"""

        for idx, line in enumerate(self.filehandle):
            currentLine = line.strip()
            if not currentLine:
                pass
            elif re.match('^#', currentLine):
                pass
            else:
		m = re.search(r"^(.*)\t(.*)\t(.*)\t(.*)$", currentLine)
		if m == None:
			raise Exception("Error bad format for line:{} in {}".format(idx,self.inputFile))

                if m.group(1) not in self.dConditions:
                    self.dConditions[m.group(1)] = []
                self.dConditions[m.group(1)].append((m.group(2),m.group(3),m.group(4)))

        return self.dConditions

    def getlLabels(self):
        """return list of labels"""

        lLabels = []
        for cond in sorted(self.dConditions.iterkeys()):
            for entry in self.dConditions[cond]:
                lLabels.append(entry[0])
        return lLabels 

    def getlConditions(self):
        """return entries ordered per condition"""

        lConditions = []
        idx = 0
        for i,cond in enumerate(sorted(self.dConditions.iterkeys())):
            lConditions.append([])
            for entry in self.dConditions[cond]:
                lConditions[i].append(idx)
                idx += 1
        return lConditions

    def getlNucFiles(self):
        """return list of nucleosome files"""

        lNucFiles = []
        for cond in sorted(self.dConditions.iterkeys()):
            for entry in self.dConditions[cond]:
                lNucFiles.append(entry[1])
        return lNucFiles 

    def getlWigFiles(self):
        """return list of wig files"""
        
        lWigFiles = []
        for cond in sorted(self.dConditions.iterkeys()):
            for entry in self.dConditions[cond]:
                lWigFiles.append(entry[2])
        return lWigFiles 

