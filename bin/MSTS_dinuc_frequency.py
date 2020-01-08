#!/usr/bin/env python

import logging
import sys
import argparse
import math
import numpy as np
from scipy import stats
import CMSTS

from pyfaidx import Fasta
import pyBigWig

from MSTS.version import __version__
from MSTS.Graphics import Graphics

from extlib.detect_peaks import detect_peaks

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
    #print len(lRhs)
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
    parser.add_argument("--pAutocor",help="print AT and GC correlograms", action="store_true", default=False)
    parser.add_argument("--pAutocorMix",help="print AT and GC autocorrelations on a single correlogram", action="store_true", default=False)
    parser.add_argument("-ami","--autocorMin", help="start for autocorrelation analysis, default=[5]", type=int, default=5)
    parser.add_argument("-amx","--autocorMax", help="stop for autocorrelation analysis, default=[35]", type=int, default=35)
    parser.add_argument("--regression",help="detect peaks and perform a regression on autocorrelation curves. Regression curve drawn on the graph", action="store_true", default=False)
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
#        if nb == 500000:
#            break
        logging.info('Retrieving sequences in: {}'.format(seq))
        if seq not in bb.chroms().keys():
            continue
        lEntries = bb.entries(seq, 0, bb.chroms(seq))
        for entry in lEntries:
  #      for entry in lEntries[0:200000]:
          #  print entry[0] ,entry[1]
            #middle = entry[0]+(entry[1]-1-entry[0])/2
            #bug
            middle = entry[0]+math.floor((entry[1]-entry[0])/2)
            start = middle - args.distance
            end = middle + args.distance + 1
         #   print start, end, middle
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

#            if nb == 500000:
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
    ATmean = np.mean(lFreqATs[1:-1])
    GCmean = np.mean(lFreqGCs[1:-1])

    lFreqATNormalized = [x/ATmean for x in lFreqATs]
    lFreqGCNormalized = [x/GCmean for x in lFreqGCs]

    lFreqATs[0] = "NaN"
    lFreqGCs[0] = "NaN"
    lFreqATNormalized[0] = "NaN"
    lFreqGCNormalized[0] = "NaN"
    lFreqATs[-1] = "NaN"
    lFreqGCs[-1] = "NaN"
    lFreqATNormalized[-1] = "NaN"
    lFreqGCNormalized[-1] = "NaN"

    logging.info('{} entries analyzed'.format(nb))

    if args.flush:
        print("pos\tATfreq\tGCfreq\tATfreqNorm\tGCfreqNorm")
        for i in range(-args.distance,args.distance+1):
            print("{}\t{}\t{}\t{}\t{}".format(i,lFreqATs[i+args.distance],lFreqGCs[i+args.distance],lFreqATNormalized[i+args.distance],lFreqGCNormalized[i+args.distance]))


    lGrid = []
    for i in range(-args.distance+1,args.distance):
        if not i % 10:
            lGrid.append(i)

    if args.pFreq:
        Graphics.plotDistribution(range(-args.distance+1,args.distance),lFreqATs[1:-1],out="{}_AT.png".format(args.prefix), title="Frequency of dinucleotides AA,AT,TA,TT", xax="position in bp", yax="frequency", legend=["AA/AT/TA/TT"],grid=lGrid)
        Graphics.plotDistribution(range(-args.distance+1,args.distance),lFreqGCs[1:-1],out="{}_GC.png".format(args.prefix), title="Frequency of dinucleotides GG,GC,CG,CC", xax="position in bp", yax="frequency", color='green', legend=["GG/GC/CG/CC"],grid=lGrid)
    if args.pFreqNorm:
        Graphics.plotDistribution(range(-args.distance+1,args.distance),lFreqATNormalized[1:-1],out="{}_AT_Normalized.png".format(args.prefix),title="Normalized frequency of dinucleotides AA,AT,TA,TT", xax="position in bp", yax="Normalized frequency", legend=["AA/AT/TA/TT"],grid=lGrid)
        Graphics.plotDistribution(range(-args.distance+1,args.distance),lFreqGCNormalized[1:-1],out="{}_GC_Normalized.png".format(args.prefix),title="Normalized frequency of dinucleotides GG,GC,CG,CC", xax="position in bp", yax="Normalized frequency", color='green', legend=["GG/GC/CG/CC"],grid=lGrid)
    if args.pFreqNormMix:
        Graphics.plotMultiDistribution(range(-args.distance+1,args.distance),[lFreqATNormalized[1:-1],lFreqGCNormalized[1:-1]],out="{}_ATGC_Normalized.png".format(args.prefix),title="Normalized frequency of dinucleotides", xax="position in bp", yax="Normalized frequency",legend=["AA/AT/TA/TT","GG/GC/CG/CC"],color=['blue','green'], grid=lGrid)
    if args.pAutocor or args.pAutocorMix:
        lRhAT = autocorrelation(range(0,args.autocorMax-args.autocorMin),lFreqATs[args.autocorMin+args.distance:args.autocorMax+args.distance])
        lRhGC = autocorrelation(range(0,args.autocorMax-args.autocorMin),lFreqGCs[args.autocorMin+args.distance:args.autocorMax+args.distance])

        if args.pAutocor:
            Graphics.plotDistribution(range(args.autocorMin,args.autocorMax),lRhAT,out="{}_AT_Correlogram.png".format(args.prefix),title="autocorrelation of dinucleotides AA,AT,TA,TT", xax="position in bp", yax="Rh, autocorrelation coefficient", legend=["AA/AT/TA/TT"])
            Graphics.plotDistribution(range(args.autocorMin,args.autocorMax),lRhGC,out="{}_GC_Correlogram.png".format(args.prefix),title="autocorrelation of dinucleotides GG,GC,CG,CC", xax="position in bp", yax="Rh, autocorrelation coefficient", color='green', legend=["GG/GC/CG/CC"])
        if args.pAutocor and args.regression:

            logging.info("Performing peak detection")
            # detect all peaks and plot data
            lPeaks = detect_peaks(lRhAT, mpd=4)
            if (len(lPeaks) > 1):
                logging.info("Multi peaks detected - performing regression on AT")
                slope, intercept, r_value, p_value, std_err = stats.linregress([x for x,y in enumerate(lPeaks)],lPeaks)

                lPhasePeriods = []
                for i in range(1,len(lPeaks)):
                    lPhasePeriods.append(lPeaks[i]-lPeaks[i-1])
                    meanlPhasePeriods = np.mean(lPhasePeriods)
                    stdlPhasePeriods = np.std(lPhasePeriods)

                logging.info("Equation: x*{} + {}".format(slope,intercept))
                logging.info("R2: {}".format(r_value))
                logging.info("p-value: {}".format(p_value))

            Graphics.plotDistributionWithRegressionDinuc(range(args.autocorMin,args.autocorMax),lRhAT,[x for x,y in enumerate(lPeaks)],lPeaks,slope, intercept,r_value,p_value,meanlPhasePeriods,stdlPhasePeriods,out="{}_AT_Correlogram.png".format(args.prefix),title="autocorrelation of dinucleotides AA,AT,TA,TT", xax="position in bp", yax="Rh, autocorrelation coefficient")

            lPeaks = detect_peaks(lRhGC, mpd=4)
            if (len(lPeaks) > 1):
                logging.info("Multi peaks detected - performing regression on AT")
                slope, intercept, r_value, p_value, std_err = stats.linregress([x for x,y in enumerate(lPeaks)],lPeaks)

                lPhasePeriods = []
                for i in range(1,len(lPeaks)):
                    lPhasePeriods.append(lPeaks[i]-lPeaks[i-1])
                    meanlPhasePeriods = np.mean(lPhasePeriods)
                    stdlPhasePeriods = np.std(lPhasePeriods)

                logging.info("Equation: x*{} + {}".format(slope,intercept))
                logging.info("R2: {}".format(r_value))
                logging.info("p-value: {}".format(p_value))

            Graphics.plotDistributionWithRegressionDinuc(range(args.autocorMin,args.autocorMax),lRhGC,out="{}_GC_Correlogram.png".format(args.prefix),title="autocorrelation of dinucleotides GG,GC,CG,CC", xax="position in bp", yax="Rh, autocorrelation coefficient", color='green')

        if args.pAutocorMix:
            Graphics.plotMultiDistribution(range(args.autocorMin,args.autocorMax),[lRhAT,lRhGC],out="{}_ATGC_Correlogram.png".format(args.prefix),title="autocorrelation of dinucleotides", xax="position in bp", yax="Rh, autocorrelation coefficient",legend=["AA/AT/TA/TT","GG/GC/CG/CC"],color=['blue','green'])
