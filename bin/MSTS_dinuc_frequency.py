#!/usr/bin/env python

import logging
import sys
import argparse

import numpy as np
import CMSTS

from pyfaidx import Fasta
import pyBigWig

from MSTS.version import __version__
from MSTS.Graphics import Graphics

def autocorrelation(lXs,lYs):
    lRhs = []
    Ym = np.mean(lYs)
    N = len(lYs)
    C0 = 0.0
    for i in lYs:
        C0 += (i-Ym)*(i-Ym)
    C0 = C0/N
    for i in range(0,N):
        Ch = 0.0
        for j in range(0,(N-i-1)):
            Ch += (lYs[j]-Ym)*(lYs[j+i]-Ym)
        Ch = Ch/N
        lRhs.append(Ch/C0)
    print len(lRhs)
    return lRhs




if __name__ == '__main__':

    program = sys.argv[0]
    version = __version__
    description = 'TODO'

    "les fragments sont consideres bien positionnes, si ils ne sont pas assez long, la sequence est allongee de chaque cote pour respecter la distance demandee. IL faut donc des fragments centres. si ce sont des reads simples, il faut les etendre de faon a considerer la sequence voulue... jusqu a 147 par exemple."

    parser = argparse.ArgumentParser(prog=program)
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("--version", action='version', version='{} {}'.format(program,version))
    parser.add_argument("FastaFile", help="Genome in fasta format", type=str)
    parser.add_argument("BigBedFile", help="nucleosome positions/centers in bigBed file ", type=str)
    parser.add_argument("-d","--distance", help="distance to the center", type=int, default=70)
    parser.add_argument("--flush",help="print frequencies on stdout to save to a file, > frequencies.out", action="store_true", default=False)
    parser.add_argument("--pFreq",help="print AT and GC frequencies plots", action="store_true", default=False)
    parser.add_argument("-p","--prefix", help="prefix use for plot file output name", type=str, default="freq")
    parser.add_argument("--pFreqNorm",help="print AT and GC Normalized frequencies plots", action="store_true", default=False)
    parser.add_argument("--pFreqNormMix",help="print AT and GC Normalized frequencies on single plot", action="store_true", default=False)
    parser.add_argument("--pAutocor",help="print AT and GC autocorrelation plots", action="store_true", default=False)
    parser.add_argument("--pAutocorMix",help="print AT and GC autocorrelation on single plot", action="store_true", default=False)
    parser.add_argument("-ami","--autocorMin", help="start for autocorrelation analysis", type=int, default=5)
    parser.add_argument("-amx","--autocorMax", help="stop for autocorrelation analysis", type=int, default=35)
    parser.add_argument("-b","--buffer", help="size of chunk (nb sequences) to keep in memory before analysis", type=int, default=1000000)
    parser.add_argument("-v", "--verbosity", type=int, choices=[1,2,3],
                        help="increase output verbosity 1=error, 2=info, 3=debug")

    args = parser.parse_args()

    logLevel='ERROR'
    if args.verbosity == 1:
        logLevel = 'ERROR'
    if args.verbosity == 2:
        logLevel = 'INFO'
    if args.verbosity == 3:
        logLevel = 'DEBUG'
    logging.getLogger().setLevel(logLevel)

    if args.autocorMin < -args.distance:
        logging.error("autcorMin value too small: {}, minimum is: {}".args.autocorMin, -args.distance)
        sys.exit(1) 
    if args.autocorMax > args.distance:
        logging.error("autcorMax value too high: {}, maximum is: {}".args.autocorMax, args.distance)
        sys.exit(1) 


    lFreqATs = [0.0]*(args.distance*2+1)
    lFreqGCs = [0.0]*(args.distance*2+1)
    nb = 0
    nbEntries = 0
    buffSize = args.buffer
    lEntriesToAnalyze = []
    bb = pyBigWig.open(args.BigBedFile)
    sequences = Fasta(args.FastaFile)
    for seq in sequences.keys():
#        if nb == 50000:
#            break
        logging.info('Retrieving sequences in: {}'.format(seq))
        if seq not in bb.chroms().keys():
            continue
        lEntries = bb.entries(seq, 0, bb.chroms(seq))
        for entry in lEntries:
            middle = entry[0]+(entry[1]-entry[0])/2
            start = middle - args.distance
            end = middle + args.distance + 1
            if (start < 0 or end > bb.chroms(seq)):
                logging.info('One entry out of boundaries, not take into account')
                continue
            eSeq = sequences[seq][start:end]
            lEntriesToAnalyze.append(str(eSeq))
            nbEntries += 1
            nb +=1
            if nbEntries == buffSize:
                lATs, lGCs = CMSTS.dinucfrequency(args.distance,lEntriesToAnalyze)
                nbEntries = 0
                lEntriesToAnalyze = []
                for i,val in enumerate(lATs):
                   lFreqATs[i] += lATs[i]
                   lFreqGCs[i] += lGCs[i]

            if not nb%buffSize:
                logging.info("Computing frequencies, nb sequences: {}".format(nb))
               
#            if nb == 50000:
#                break

    if nbEntries > 0:
        lATs,lGCs = CMSTS.dinucfrequency(args.distance,lEntriesToAnalyze)
        nbEntries = 0
        lEntriesToAnalyze = []
        for i,val in enumerate(lATs):
            lFreqATs[i] += lATs[i]
            lFreqGCs[i] += lGCs[i]

    lFreqATs = [ x/nb for x in lFreqATs ]
    lFreqGCs = [ x/nb for x in lFreqGCs ]
    ATmean = np.mean(lFreqATs)
    GCmean = np.mean(lFreqGCs)

    lFreqATNormalized = [x/ATmean for x in lFreqATs]
    lFreqGCNormalized = [x/GCmean for x in lFreqGCs]

    logging.info('{} positions analyzed'.format(nb))

    if args.flush:
        print "pos\tATfreq\tGCfreq\ATfreqNorm\tGCfreqNorm"
        for i in range(0,2*args.distance+1):
            print "{}\t{}\t{}\t{}".format(lFreqATs[i],lFreqGCs[i],lFreqATNormalized,lFreqGCNormalized)

    if args.pFreq:
        Graphics.plotDistribution(range(-args.distance+1,args.distance),lFreqATs[1:-1],out="{}AT.png".format(args.prefix), title="Frequency of dinucleotides AA,AT,TA,TT", xax="position in bp", yax="frequency", legend=["AA/AT/TA/TT"])
        Graphics.plotDistribution(range(-args.distance+1,args.distance),lFreqGCs[1:-1],out="{}GC.png".format(args.prefix), title="Frequency of dinucleotides GG,GC,CG,CC", xax="position in bp", yax="frequency", color='green', legend=["GG/GC/CG/CC"])
    if args.pFreqNorm:
        Graphics.plotDistribution(range(-args.distance+1,args.distance),lFreqATNormalized[1:-1],out="{}ATNormalized.png".format(args.prefix),title="Normalized frequency of dinucleotides AA,AT,TA,TT", xax="position in bp", yax="Normalized frequency", legend=["AA/AT/TA/TT"])
        Graphics.plotDistribution(range(-args.distance+1,args.distance),lFreqGCNormalized[1:-1],out="{}GCNormalized.png".format(args.prefix),title="Normalized frequency of dinucleotides GG,GC,CG,CC", xax="position in bp", yax="Normalized frequency", color='green', legend=["GG/GC/CG/CC"])
    if args.pFreqNormMix:
        Graphics.plotMultiDistribution(range(-args.distance+1,args.distance),[lFreqATNormalized[1:-1],lFreqGCNormalized[1:-1]],out="{}ATGCNormalized.png".format(args.prefix),title="Normalized frequency of dinucleotides", xax="position in bp", yax="Normalized frequency",legend=["AA/AT/TA/TT","GG/GC/CG/CC"],color=['blue','green'])
    if args.pAutocor or args.pAutocorMix:
        lRhAT = autocorrelation(range(0,args.autocorMax-args.autocorMin),lFreqATs[args.autocorMin+args.distance:args.autocorMax+args.distance])
        lRhGC = autocorrelation(range(0,args.autocorMax-args.autocorMin),lFreqGCs[args.autocorMin+args.distance:args.autocorMax+args.distance])
        if args.pAutocor: 
            Graphics.plotDistribution(range(args.autocorMin,args.autocorMax),lRhAT,out="{}ATAutocorrelation.png".format(args.prefix))
            Graphics.plotDistribution(range(args.autocorMin,args.autocorMax),lRhGC,out="{}GCAutocorrelation.png".format(args.prefix))
        if args.pAutocorMix:
            Graphics.plotMultiDistribution(range(args.autocorMin,args.autocorMax),[lRhAT,lRhGC],out="{}ATGCAutocorrelation.png".format(args.prefix),legend=["AA/AT/TA/TT","GG/GC/CG/CC"],color=['blue','green'])
        pass  
