#/usr/bin/env python

""" BamConverter: from BAM to wig/bed/cov/size format

    Your Bam file must be sorted and indexed

----------------------------------------------------------------------------------------
| export | mode            | single read        | paired reads
----------------------------------------------------------------------------------------
| bed    | single          | start to end read  | start to end for both reads
|        | fragment        | start to end read  | start read1 to end read2
|        | fragment-middle | middle of the read | middle of start read 1 to end read 2 
----------------------------------------------------------------------------------------
| cov    | single          | start to end read  | start to end for both reads
|        | fragment        | start to end read  | start read1 to end read2
|        | fragment-middle | middle of the read | middle of start read 1 to end read 2 
----------------------------------------------------------------------------------------
| wig    | single          | start to end read  | start to end for both reads
|        | fragment        | start to end read  | start read1 to end read2
|        | fragment-middle | middle of the read | middle of start read 1 to end read 2 
----------------------------------------------------------------------------------------
| size   | single          | start to end read  | start to end for both reads
|        | fragment        | start to end read  | start read1 to end read2
|        | fragment-middle | middle of the read | start read 1 to end read 2 
---------------------------------------------------------------------------------------

Note: the middle of a read/fragment is one base for odd read length 
      and two bases for even pair read length

"""

import sys
import os
#import time
import math
import logging
import pysam
import pyBigWig
import numpy as np
import array

class BamConverter(object):

    def __init__(self, BamFile='', mode='single', prefix='out', bed=False, wig=False, cov=False, size=False, genome='', window=0, keepPosBigBedFile='', minFragSize=0, maxFragSize=2000, minDepCov=0, maxDepCov=1000000, logLevel='ERROR'):
        """constructor"""

        self.BamFile = BamFile
        self.mode = mode
        self.prefix = prefix
        self.bed = bed
        self.wig = wig
        self.cov = cov
        self.size = size
        self.genome = genome
        self.window = window
        self.keepPosBigBedFile = keepPosBigBedFile
        self.minFragSize = minFragSize
        self.maxFragSize = maxFragSize
        self.minDepCov = minDepCov
        self.maxDepCov = maxDepCov
        self.logLevel = logLevel
        logging.basicConfig(level=self.logLevel)

        try:
            self.bamFileH = pysam.AlignmentFile(self.BamFile, "r", check_sq=True)
        except:
            logging.error("Can not read input file: {}".format(self.BamFile))

        self.dReferences = {k: v for k, v in zip(self.bamFileH.references, self.bamFileH.lengths)}

        self.lRefOrder = [x for x in self.dReferences] 

        if self.genome:
            self.lRefOrder = self.getReferenceOrder()

        if self.keepPosBigBedFile:
            self.bb = pyBigWig.open(self.keepPosBigBedFile)

        if self.bed or self.wig or self.cov or self.size:
            self.checkIfOutputFileExist()

        if self.wig:
            self.printWigHeader()

        self.nbTotalFrags = 0
        self.nbAnalyzedFrags = 0


    def __del__(self):

        self.bamFileH.close()

    def convert(self):
        """convert Bam to wig/bed/cov/size"""

        dAllFragmentSize = {k: v for k, v in zip(range(1,1001), [0]*1000)}

        for seq in self.lRefOrder:
            logging.info("Processing ref: {}".format(seq))

            currentSeq = array.array('i')
            for i in range(0,self.dReferences[seq]):
                currentSeq.append(0)

            if self.mode == "single" or self.mode == "single-expanded":
                currentSeq, dFragmentSize, lBedTracks = self.convertWithSingleMode(seq, currentSeq)
            elif self.mode == "fragment" or self.mode == "fragment-middle":
                currentSeq, dFragmentSize, lBedTracks = self.convertWithFragmentMode(seq, currentSeq)

            if self.wig:
                self.printWig(seq,currentSeq)
            if self.cov: 
                self.printCov(seq,currentSeq)
            if self.bed:
                self.printBed(lBedTracks)

            for size in dAllFragmentSize:
                if size in dFragmentSize:
                    dAllFragmentSize[size] = dAllFragmentSize[size] + dFragmentSize[size]

        if self.size:
            self.printFragmentSize(dAllFragmentSize)
        logging.info("Total Nb Fragments: {} -- Exported Nb Fragments: {}".format(self.nbTotalFrags, self.nbAnalyzedFrags))

    def isFragmentFullyIncludeInABedTrack(self, seq, start, end):
        "..."
        try:
            lEntries = self.bb.entries(seq, start, end)
            if lEntries:
                if lEntries[0][0] <= start and lEntries[0][1] >= end:
                    return True
            else:
                return False
        except Exception as e:
            logging.debug('cannot retrieve: {}-{}-{}, message:{}'.format(seq, start, end, e))
            return False

    def convertWithSingleMode(self, seq, currentSeq):
        """convert in single mode"""
        
        lBedTracks = []
        dFragmentSize = {}
        nbTotalFrags = 0
        nbAnalyzedFrags = 0

        for read in self.bamFileH.fetch(reference=seq):
            # process single and paired read in same manner
            nbTotalFrags+=1 

            start = read.get_reference_positions()[0]
            end = read.get_reference_positions()[-1]
            
            if self.mode == "single-expanded":
                if read.is_reverse:
                    start = max(0,end - 146) 
                else:
                    end = min(len(currentSeq)-1,start + 146)
            # to prevent bug with mapping close to start/end 
            else:
                start = max(0,start) 
                end = min(len(currentSeq)-1,end)
 
                

            if start and end:
                if self.keepPosBigBedFile:
                    if not self.isFragmentFullyIncludeInABedTrack(seq,start,end):
                        continue

                if (end-start+1) < self.minFragSize or (end-start+1) > self.maxFragSize:
                    continue

                if end-start+1 in dFragmentSize:
                    dFragmentSize[end-start+1] += 1
                else:
                    dFragmentSize[end-start+1] = 1
 
                nbAnalyzedFrags+=1

                lBedTracks.append((seq,start,end))

                for x in range(start,end+1):
                    #print end+1
                    #print len(currentSeq)
                    #print read.reference_start
                    #print read.reference_end
                    #print read.query_name
                    currentSeq[x] += 1
            
        if self.minDepCov != 0 or self.maxDepCov != 1000000:
            currentSeq, lBedTracks, dFragmentSize = self.reduceDataToSpecificDepCov(currentSeq, lBedTracks, dFragmentSize)

        logging.info("{}: Total Nb Fragments: {} -- Exported Nb Fragments: {}".format(seq, nbTotalFrags, nbAnalyzedFrags))
        self.nbTotalFrags += nbTotalFrags
        self.nbAnalyzedFrags += nbAnalyzedFrags

        return currentSeq, dFragmentSize, lBedTracks 


    def convertWithFragmentMode(self, seq, currentSeq):
        """convert in fragment/fragment-middle mode"""

        lBedTracks = []
        dFragmentSize = {}
        nbTotalFrags = 0
        nbAnalyzedFrags = 0

        for read in self.bamFileH.fetch(reference=seq):
 
            start = None
            end = None
            # process paired reads
            if read.is_paired:
                if read.is_read1:
                    nbTotalFrags+=1 
                    if (read.template_length < 0):
                        start = read.next_reference_start
                        end = read.next_reference_start-read.template_length-1
                    else:
                        start = read.reference_start
                        end = read.reference_start + read.template_length-1 
         #           nbr += 1
            # process single
            else:
                start = read.get_reference_positions()[0]
                end = read.get_reference_positions()[-1]
                    
            if start and end:
                if self.keepPosBigBedFile:
                    if not self.isFragmentFullyIncludeInABedTrack(seq,start,end):
                        continue

                if (end-start+1) < self.minFragSize or (end-start+1) > self.maxFragSize:
                    continue

                if end-start+1 in dFragmentSize:
                    dFragmentSize[end-start+1] += 1
                else:
                    dFragmentSize[end-start+1] = 1
                ##debug
         #       print "{}\t{}\t{}".format(read.query_name,start,end)


                ##debug 
                nbAnalyzedFrags+=1 
                if self.mode == "fragment":
                    lBedTracks.append((seq,start,end))
                    for x in range(start,end+1):
                        currentSeq[x] += 1
                elif self.mode =="fragment-middle":
                    #middle = end-((end-start+1)/2)
                    middle = int(start+math.ceil((end-start+1)/2.0))
                    ##debug
         #           if middle == 198:
         #               print "198: {}".format(read.query_name)
         #           if middle == 199:
         #               print "199: {}".format(read.query_name)
                    ##debug
                    window = self.window

                    if (end-start+1) < self.window:
                        if (end-start+1)%2:
                            window = (end-start+1)/2-1
                        else:
                            window = (end-start+1)/2
                    
                 #   if (end-start+1)%2:
                        
                 #       for x in range(middle-window, middle+1+window):
                 #           currentSeq[x] += 1 
                 #           nbmiddle += 1
                 #       lBedTracks.append((seq,middle - window ,(middle+1) + window))
                #    else:
                 #       for x in range(middle-window,middle+window):
                 #           currentSeq[x] += 1
                 #           nbmiddle +=1
                 #       lBedTracks.append((seq,middle - window,middle + window))
                    
        #            nbmiddle = 0
                    for x in range(middle-window, middle+1+window):
                        try:
                            currentSeq[x] += 1 
                        except:
                            print(x)
        #                nbmiddle += 1
        #            print nbmiddle
                    lBedTracks.append((seq,middle - window , middle+1+ window))

        if self.minDepCov != 0 or self.maxDepCov != 1000000:
            currentSeq, lBedTracks, dFragmentSize = self.reduceDataToSpecificDepCov(currentSeq, lBedTracks, dFragmentSize)

        logging.info("{}: Total Nb Fragments: {} -- Exported Nb Fragments: {}".format(seq, nbTotalFrags, nbAnalyzedFrags))
        self.nbTotalFrags += nbTotalFrags
        self.nbAnalyzedFrags += nbAnalyzedFrags

        return currentSeq, dFragmentSize, lBedTracks 


    def reduceDataToSpecificDepCov(self, lCov, lBedTracks, dFragmentSize):
        """reduce data to min/max depth of coverage"""

        logging.info("Filtering results with depth of coverage: {} < x < {}".format(self.minDepCov, self.maxDepCov))
        lToFilteredOut = [1 if x < self.minDepCov or x > self.maxDepCov else 0 for x in lCov]
        lCovFilteredOut  = ['F' if x < self.minDepCov or x > self.maxDepCov else x for x in lCov]
        lBedTracksFilteredOut = []
        for i,fragment in enumerate(lBedTracks):
            if sum(lToFilteredOut[fragment[1]:fragment[2]+1]) == 0:
                lBedTracksFilteredOut.append(lBedTracks[i])
            else:
                dFragmentSize[fragment[2]-fragment[1]+1] -= 1
                
        return lCovFilteredOut, lBedTracksFilteredOut, dFragmentSize 


    def getReferenceOrder(self):
        """get reference order from file"""

        lReferences = []
        logging.info("Reading genome file: {}".format(self.genome))   
        try:
            with open(self.genome, 'r') as gFile:
                for line in gFile:
                    values = line.rstrip().split("\t")
                    if len(values) != 2:
                        logging.error("Genome File not well formatted. \
                                       Expected format: <refname><TAB><size>")
                    else:
                        lReferences.append(values[0])
        except:
            logging.error("Can not read/open genome file: {}".format(self.genome))
            sys.exit(1)

        return lReferences

    def checkIfOutputFileExist(self):
        """check if required outputs exist"""

        if self.size:
            self._checkAndRemoveFile(".size")
        if self.cov:
            self._checkAndRemoveFile(".cov")
        if self.bed:
            self._checkAndRemoveFile(".bed")
        if self.wig:
            self._checkAndRemoveFile(".wig")

    def _checkAndRemoveFile(self, suffix):
        """check if file exists and remove it"""

        if os.path.isfile(self.prefix + suffix):
            logging.info("Removing existing file: {}{}".format(self.prefix,suffix))
            os.remove(self.prefix + suffix)


    def printFragmentSize(self, dFragmentSize):
        """write output in size format"""
        
        try:
            with open(self.prefix + ".size", 'a') as f:
                f.write("size\tnb fragments\n")
                for size,value in dFragmentSize.items():
                    f.write("{}\t{}\n".format(size, value))
            f.close()
        except:
            logging.error("Can not write results in {}.size".format(self.prefix))


    def printBed(self, lBedTracks):
        """write output in bed format"""

        try:
            with open(self.prefix + ".bed", 'a') as f:
                for fragment in lBedTracks:
                    #f.write("{}\t{}\t{}\n".format(fragment[0],fragment[1],fragment[2]+1))
                    f.write("{}\t{}\t{}\n".format(fragment[0],fragment[1],fragment[2]))
            f.close()
        except:
            logging.error("Can not write results in {}.bed".format(self.prefix))


    def printCov(self, seq, lcov):
        """write output in bed format"""

        try:
            with open(self.prefix + ".cov", 'a') as f:
                for pos,cov in enumerate(lcov):
                    f.write("{}\t{}\t{}\n".format(seq,pos+1,cov))
            f.close()
        except:
            logging.error("Can not write results in {}.cov".format(self.prefix))


    def printWigHeader(self):
        """write wig header"""

        try:
            with open(self.prefix + ".wig", 'a') as f:
                f.write("track type=wiggle_0 name=fileName description=fileName\n")
            f.close()
        except:
            logging.error("Can not write results in {}.wig".format(self.prefix))

    def printWig(self, seq, lcov):
        """write output in wig format"""

        try:
            with open(self.prefix + ".wig", 'a') as f:
                f.write("variableStep chrom={}\n".format(seq))
                if self.keepPosBigBedFile:
                    for pos,cov in enumerate(lcov):
                        if self.isFragmentFullyIncludeInABedTrack(seq,pos,pos+1):
                            f.write("{}\t{}\n".format(pos+1,cov))
                else:
                    for pos,cov in enumerate(lcov):
                        f.write("{}\t{}\n".format(pos+1,cov))
            f.close()
        except:
            logging.error("Can not write results in {}.wig".format(self.prefix))


