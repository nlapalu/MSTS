#!/usr/bin/env python

import sys
import logging
import argparse
import pyBigWig
import math

import CMSTS
import numpy as np
from scipy import stats

from MSTS.version import __version__
from MSTS.Graphics import Graphics

from extlib.detect_peaks import detect_peaks

if __name__ == "__main__":

    program = sys.argv[0]
    version = __version__
    description = 'Draw phasogram from bigWig'

    parser = argparse.ArgumentParser(prog=program)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--version', action='version', version='{} {}'.format(program,version))

    parser.add_argument("bigWig", help="Input bigWig File", type=str)
    parser.add_argument("-w", "--window", help="window size to compute phases", type=int, default=1200)
    parser.add_argument("--flush",help="print phases on stdout to save in file, > phases.out", action="store_true", default=False)
    parser.add_argument("--regression",help="detect peaks and perform a regression. Regression curve drawn on the graph", action="store_true", default=False)
    parser.add_argument("--norm", help="normalize the signal with the mean coverage", action="store_true", default=False)
    parser.add_argument("-o", "--out", help="name of output graph", type=str, default="graph.png")
    parser.add_argument("-t", "--title", help="title text", type=str, default=None)
    parser.add_argument("-x", "--xax", help="x axis text", type=str, default="window, bp")
    parser.add_argument("-y", "--yax", help="y axis text", type=str, default="signal coverage")
    parser.add_argument("-b", "--bigBed", help="bigBedFile, use to limit phasogram to specific regions", type=str, default=None)

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

    lPhases = [0]*args.window
    buffSize = 1000000
    debug = 0
    nb_windows = 0
    allvalues = []

    logging.info('Buffer size: {} bases'.format(buffSize))

    bw = pyBigWig.open(args.bigWig)
    for chrom in bw.chroms():
        logging.info('Reading sequence: {}'.format(chrom))
        if args.bigBed:
            bb = pyBigWig.open(args.bigBed)
            if chrom not in bb.chroms().keys():
                continue
            lEntries = bb.entries(chrom, 0, bw.chroms(chrom))
            for entry in lEntries:
                start = entry[0]
                loop = 0
                while (start < entry[1]):
                    stop = min(start+buffSize, entry[1])
                    logging.info('Requesting values: {}:{}-{}'.format(chrom,start,stop))
                    values = bw.values(chrom, start, stop)
                    iStart = 0
                    iEnd = len(values)
                    for idx,val in enumerate(values):
                        #if math.isnan(val) and iStart != None:
                        if (math.isnan(val) or val==0) and iStart != None:
                            iEnd = idx
                            if len(values[iStart:iEnd]) > args.window:
                                l = CMSTS.phasogram(args.window,[float(i) for i in values[iStart:iEnd]])
                                if args.norm:
                                    nb_windows += len(values[iStart:iEnd]) - args.window + 1
                                    allvalues.extend(values[iStart:iEnd])
                                lPhases = [a1 + b1 for a1, b1 in zip(l,lPhases)]
                            iStart = None
                        elif (not math.isnan(val) or val!=0) and iStart == None:
                            iStart = idx
                        else:
                            iEnd = idx
                    if iStart != None:
                        if len(values[iStart:iEnd]) > args.window:
                            l = CMSTS.phasogram(args.window,[float(i) for i in values[iStart:iEnd]])
                            if args.norm:
                                nb_windows += len(values[iStart:iEnd]) - args.window + 1
                                allvalues.extend(values[iStart:iEnd])
                            lPhases = [a1 + b1 for a1, b1 in zip(l,lPhases)]
                    loop+=1
                    start = stop
        else:
            start = 0
            loop = 0
            while (start < bw.chroms(chrom)):
                stop = min(start+buffSize, bw.chroms(chrom))
                logging.info('Requesting values: {}:{}-{}'.format(chrom,start,stop))
                values = bw.values(chrom, start, stop)
                iStart = 0
                iEnd = len(values)

                for idx,val in enumerate(values):
                    if math.isnan(val) and iStart != None :
                        iEnd = idx
                        if len(values[iStart:iEnd]) > args.window:
                            l = CMSTS.phasogram(args.window,[float(i) for i in values[iStart:iEnd]])
                            if args.norm:
                                nb_windows += len(values[iStart:iEnd]) - args.window + 1
                                allvalues.extend(values[iStart:iEnd])
                            lPhases = [a1 + b1 for a1, b1 in zip(l,lPhases)]
                        iStart = None
                    elif not math.isnan(val) and iStart == None:
                        iStart = idx
                    else:
                        iEnd = idx
                if iStart != None:
                    if len(values[iStart:iEnd]) > args.window:
                        l = CMSTS.phasogram(args.window,[float(i) for i in values[iStart:iEnd]])
                        if args.norm:
                            nb_windows += len(values[iStart:iEnd]) - args.window + 1
                            allvalues.extend(values[iStart:iEnd])
                        lPhases = [a1 + b1 for a1, b1 in zip(l,lPhases)]
                loop+=1
                start = stop

    if args.norm:
        mean = sum(allvalues)/len(allvalues)
        lPhases = [v/(nb_windows*mean) for v in lPhases]

    title = "phasogram of {}".format(args.bigWig)
    if args.title:
        title = args.title

    if args.regression:
        logging.info("Performing peak detection")
        # detect all peaks and plot data
        lPeaks = detect_peaks(lPhases, mpd=147)
        if (len(lPeaks) > 1):
            logging.info("Multi peaks detected - performing regression")
            slope, intercept, r_value, p_value, std_err = stats.linregress([x for x,y in enumerate(lPeaks)],lPeaks)

            lPhasePeriods = []
            for i in range(1,len(lPeaks)):
                lPhasePeriods.append(lPeaks[i]-lPeaks[i-1])
                meanlPhasePeriods = np.mean(lPhasePeriods)
                stdlPhasePeriods = np.std(lPhasePeriods)

            logging.info("Equation: x*{} + {}".format(slope,intercept))
            logging.info("R2: {}".format(r_value))
            logging.info("p-value: {}".format(p_value))

            logging.info("Drawing graph in {}".format(args.out))
            Graphics.plotDistributionWithRegression([x for x in range(0,args.window)],lPhases, [x for x,y in enumerate(lPeaks)],lPeaks,slope, intercept,r_value,p_value,meanlPhasePeriods,stdlPhasePeriods, out=args.out, title=title, xax=args.xax, yax=args.yax)
        elif (len(lPeaks) == 1):
            logging.info("Only one peak detected - no regression")
            logging.info("Drawing graph in {}".format(args.out))
            Graphics.plotDistributionWithRegression([x for x in range(0,args.window)],lPhases, [x for x,y in enumerate(lPeaks)],lPeaks,0.0, 0.0,0.0,0.0,lPeaks[0],0, out=args.out, title=title, xax=args.xax, yax=args.yax)
        else:
            logging.info("No peak detected - drawing simple graph")
            logging.info("Drawing graph in {}".format(args.out))
            Graphics.plotDistribution([x for x in range(0,args.window)],lPhases,out=args.out, title=title, xax=args.xax, yax=args.yax)
    else:
        logging.info("Drawing graph in {}".format(args.out))
        Graphics.plotDistribution([x for x in range(0,args.window)],lPhases,out=args.out, title=title, xax=args.xax, yax=args.yax)

    if args.flush:
        for x in range(0,args.window):
            print("{}\t{}".format(x,lPhases[x]))

