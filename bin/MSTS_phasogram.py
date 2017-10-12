#/usr/bin/env python

import sys
import logging
import argparse
import pyBigWig
import math

import CMSTS

from MSTS.version import __version__
from MSTS.Graphics import Graphics

if __name__ == "__main__":

    program = sys.argv[0]
    version = __version__
    description = 'Draw phasogram from bigWig'

    parser = argparse.ArgumentParser(prog=program)
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("bigWig", help="Input bigWig File", type=str)
    parser.add_argument("-w", "--window", help="window size to compute phases", type=int, default=1200)
    parser.add_argument("--flush",help="print phases on stdout to save in file, > phases.out", action="store_true", default=False)
    parser.add_argument("-o", "--out", help="name of output graph", type=str, default="graph.png")
    parser.add_argument("-t", "--title", help="title text", type=str, default=None)
    parser.add_argument("-x", "--xax", help="x axis text", type=str, default="window, bp")
    parser.add_argument("-y", "--yax", help="y axis text", type=str, default="signal coverage")

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
    buffSize = 10000000
    debug = 0

    logging.info('Buffer size: {} bases'.format(buffSize))

    bw = pyBigWig.open(args.bigWig)
    for chrom in bw.chroms():
        logging.info('Reading sequence: {}'.format(chrom))
        start = 0
        while (start < bw.chroms(chrom)): 
            stop = min(start+buffSize - 1, bw.chroms(chrom))
            logging.info('Requesting values: {}:{}-{}'.format(chrom,start,stop))
            values = bw.values(chrom, start, stop)
            iStart = start
            iEnd = stop
          
            for idx,val in enumerate(values):
                if math.isnan(val) and iStart != None:
                    iEnd = idx-1
                    l = CMSTS.phasogram(args.window,[float(i) for i in values[iStart:iEnd]])
                    lPhases = [a1 + b1 for a1, b1 in zip(l,lPhases)]
                    iStart = None
                elif not math.isnan(val) and iStart == None:
                    iStart = idx
                else:
                    iEnd = idx
            l = CMSTS.phasogram(args.window,[float(i) for i in values[iStart:iEnd]])
            lPhases = [a1 + b1 for a1, b1 in zip(l,lPhases)]
            start = stop + 1

    title = "phasogram of {}".format(args.bigWig)
    if args.title:
        title = args.title        

    Graphics.plotDistribution([x for x in range(0,args.window)],lPhases,out=args.out, title=title, xax=args.xax, yax=args.yax)

    if args.flush:
        for x in range(0,args.window):
            print "{}\t{}".format(x,lPhases[x])
