#!/usr/bin/env python

import sys
import logging
import time
import argparse
import pyBigWig
import math

from MSTS.version import __version__
from MSTS.Parser.SimpleGffParser import SimpleGffParser
from MSTS.Db.FeatureDB import FeatureDB
from MSTS.Graphics import Graphics

def getCoordinatesFreeFromOtherGenes(lFeatures,feature,start, pivot='start'):
    "...."

    startNew = start
    endNew = start
    for feat in lFeatures:
        if feature.strand == 1:
            if pivot == 'start':
                if feat.end <= feature.start and feat.end > startNew:
                    startNew = feat.end
            elif pivot == 'end':
                if feat.start >= feature.end and feat.start < endNew:
                    endNew = feat.start
        elif feature.strand == -1:
            if pivot == 'start':
                if feat.start >= feature.end and feat.start < startNew:
                    startNew = feat.start
            elif pivot == 'end':
                if feat.end <= feature.start and feat.end > endNew:
                    endNew = feat.end

    if pivot == 'start':
        return startNew 
    elif pivot == 'end':
        return endNew


def getBasesOverlappingOtherGenes(lOverlappingFeatures,Start, End, regionStart, regionEnd):
    """...."""


    lOtherFeatures = []
    for feat in lOverlappingFeatures:
        if feat.start <= Start and feat.end >= End:
            pass
        else:
            lOtherFeatures.append(feat)

    dBases = {}
    for val in range(regionStart, regionEnd):
        dBases[val] = 0
    for feat in lOtherFeatures:
        for base in range(feat.start,feat.end):
            if base in dBases:
                dBases[base] += 1

    return [value for (key,value) in sorted(dBases.items())]


def gaussianSmoothing(data, windowWidth=3, stdev=20):
        """smoothing"""

        windowWidth = windowWidth
        stdev = stdev

        filter = [None]*(2*windowWidth*stdev+1)
        sumt = 0.0;

        for i in range(0,len(filter)):
            x = float(i - 3 * stdev)
            value =  math.exp(-(x * x) / (2 * stdev * stdev))
            filter[i] = value
            sumt += value

        for i in range(0,len(filter)):
            filter[i] /=sumt

        smoothed = [0]*len(data)
        for i in range(0,len(smoothed)-len(filter)):
            for j in range(0,len(filter)):
                smoothed[i] += data[i + j] * filter[j]

        return smoothed
        


if __name__ == "__main__":

    program = sys.argv[0]
    version = __version__
    description = 'todo, \
                   ...'

    parser = argparse.ArgumentParser(prog=program)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--version', action='version', version='{} {}'.format(program,version))

    parser.add_argument("bigWig", help="Input bigWig File", type=str)
    parser.add_argument("gff3", help="Input genome annotation in gff3 format", type=str)
    parser.add_argument("-wb","--windowBefore", help="window size to analyze before the feature, default=1000", type=int, default=1000)
    parser.add_argument("-wa","--windowAfter", help="window size to analyez after the feature, default=1000", type=int, default=1000)

    parser.add_argument("-ft","--featureType", help="feature type to analyze, default=gene", type=str, default='gene')
    parser.add_argument("-p", "--pivot", help="feature bound to use, default=start, possible values=[start,end]",type=str,default='start')
    parser.add_argument("--context", help="if set, defined features in context matter", action="store_true", default=False)
    parser.add_argument("-o", "--out", help="name of output graph", type=str, default="graph.png")
    parser.add_argument("-t", "--title", help="title text", type=str, default="title")
    parser.add_argument("-x", "--xax", help="x axis text", type=str, default="window, bp")
    parser.add_argument("-y", "--yax", help="y axis text", type=str, default="nb bases")
    parser.add_argument("-z", "--zax", help="z axis text", type=str, default="signal coverage")
    parser.add_argument("-d", "--sqliteDB", help="provide sqlite DB to avoid insertion, usefull for multi-analysis", type=str, default=None)
    parser.add_argument("-n", "--noDeleteDB", help="Do not delete SQLite DB", action="store_true", default=False)
    parser.add_argument("-s","--GaussianSmoothing", help="Perform Gaussian Smoothing on data, ", action="store_true", default=False)
    parser.add_argument("-w","--windowWidth", help="window size for Gaussian smoothing, default=3", type=int, default=3)
    parser.add_argument("-sd","--stdev", help="stdev for Gaussian smoothing, default=20", type=int, default=20)


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

    featType = args.featureType
    featTypeContext = featType
    pivot = args.pivot
    if pivot not in ['start','end']:
        logging.error("Provided pivot: '{}' not allowed, choose 'start' or 'end'".format(pivot))
        sys.exit(1)
    context = args.context

    if not args.sqliteDB:
        logging.info("Parsing gff3 gene file")
        parser = SimpleGffParser(args.gff3, logLevel)
        lFeatures = []
        for feat in parser.parse():
            lFeatures.append(feat)
        logging.info("{} features parsed".format(parser.nbFeatures))

        logging.info("Inserting features in SQLite DB")
        timestamp = int(time.time())
        db = FeatureDB('sqlite-{}.db'.format(timestamp),False,logLevel)
        db.insertlFeatures(lFeatures)
        logging.info("Insertion done")
    else:
        logging.info("Using {} file as SQLite db".format(args.sqliteDB))
        db = FeatureDB(args.sqliteDB,noCreate=True,logLevel=logLevel)

    bw = pyBigWig.open(args.bigWig)
    winBefore = args.windowBefore
    winAfter = args.windowAfter
    lPhases = [0]*(1+winBefore+winAfter)
    lPhasesNb = [0]*(1+winBefore+winAfter)
    lOtherGenesNb = [0]*(1+winBefore+winAfter)

    for chrom in db.selectReferences():
        logging.info('Requesting genes in sequence: {}'.format(chrom))
        lFeatures = db.selectFeatureTypeFromReference(chrom,featType)
        for feat in lFeatures:
            if feat.strand == 1:
                if pivot == 'start':
                    start = max(1,feat.start-winBefore)
                    end = min(bw.chroms(chrom),feat.start+winAfter)
                elif pivot == 'end':
                    start = max(1,feat.end-winBefore)
                    end = min(bw.chroms(chrom),feat.end+winAfter)
            elif feat.strand == -1:
                if pivot == 'start':
                    start = min(bw.chroms(chrom),feat.end+winBefore)
                    end = max(1,feat.end-winAfter)
                elif pivot == 'end':
                    start = min(bw.chroms(chrom),feat.start+winBefore)
                    end = max(1,feat.start-winAfter)
            else:
                logging.error("Cannot perform analysis on feature witout strand")
                sys.exit(1)

            startNew = start
            endNew = end

            if pivot == 'start':          
                lOverlappingFeatures = db.selectFeatureTypeFromCoordinates(featTypeContext,chrom,min(start,end),max(start, end))
                if context:
                    startNew = getCoordinatesFreeFromOtherGenes(lOverlappingFeatures,feat,start) 
            elif pivot == 'end':
                lOverlappingFeatures = db.selectFeatureTypeFromCoordinates(featTypeContext,chrom,min(start,end),max(start, end))
                if context:
                    endNew = getCoordinatesFreeFromOtherGenes(lOverlappingFeatures,feat,end,pivot='end') 

            index = 0
            lOtherGenesBases = []
            if feat.strand == 1:
                if pivot == 'start':
                    lValues = bw.values(chrom,startNew-1,end) 
                    decal = 0
                    if context:
                        lValues = bw.values(chrom,startNew-1,feat.end)
                        decal = (startNew-start)
                    for i in range(0+decal,min(len(lValues)+decal,winBefore+winAfter+1)):
                        lPhases[i] += lValues[i-decal]
                        lPhasesNb[i] += 1
                        index = i
                    lOtherGenesBases = getBasesOverlappingOtherGenes(lOverlappingFeatures,feat.start,feat.end,start,end)
                elif pivot == 'end':
                    
                    lValues = bw.values(chrom,start-1,endNew)
                    decal=0
                    if context:
                        lValues = bw.values(chrom,feat.start-1,endNew)
                        decal = max((winBefore-(feat.end-feat.start)),0)

                    for i in range(decal,min(decal+len(lValues),winBefore+winAfter+1)):
                        lPhases[i] += lValues[i-decal]
                        lPhasesNb[i] += 1
                        index = i
                    lOtherGenesBases = getBasesOverlappingOtherGenes(lOverlappingFeatures,feat.start,feat.end,start,end)
            elif feat.strand == -1:
                if pivot == 'start':
                    lValues = bw.values(chrom,end,startNew)
                    decal = 0
                    if context:
                        lValues = bw.values(chrom,feat.start-1,startNew)
                        decal = (startNew-start)
                    for i in range(-1+decal,max(-len(lValues)+decal-1,(-winAfter)+(-winBefore)+(-1)-1),-1):
                        lPhases[-i-1] += lValues[i-decal]
                        lPhasesNb[-i-1] += 1
                        index = i
                    lOtherGenesBases = getBasesOverlappingOtherGenes(lOverlappingFeatures,feat.start,feat.end,end,start)[::-1]
                elif pivot == 'end':
                    lValues = bw.values(chrom,endNew-1,start)
                    decal = 0
                    if context:
                        lValues = bw.values(chrom,endNew-1,feat.end)
                        decal = max((winBefore-(feat.end-feat.start)),0)
                    for i in range(-1-decal,max(-len(lValues)-decal,(-winAfter)+(-winBefore)+(-1)-1), -1):

                        lPhases[-i-1] += lValues[i+decal]
                        lPhasesNb[-i-1] += 1
                        index = i
                    lOtherGenesBases = getBasesOverlappingOtherGenes(lOverlappingFeatures,feat.start,feat.end,end,start)[::-1]
                
            else:
                pass

            for i in range(0,len(lOtherGenesBases)):
                lOtherGenesNb[i] += lOtherGenesBases[i]


    lAveragePhases = [0]*(1+winBefore+winAfter)

    for a,b in enumerate(lPhases):
        if lPhasesNb[a] != 0:
            lAveragePhases[a] = lPhases[a]/lPhasesNb[a]


    if args.GaussianSmoothing:
        logging.info("Smoothing data with Gaussian blur, window: {}, stdev: {}".format(args.windowWidth, args.stdev))
        lAveragePhases = gaussianSmoothing(lAveragePhases, args.windowWidth, args.stdev)
        lenFilter = 2*args.windowWidth*args.stdev+1

    if args.noDeleteDB: 
        logging.info('SQLite db: {} not removed'.format(db.getDbFileName()))
    else: 
        logging.info('SQLite db: {} removed'.format(db.getDbFileName()))
        db.deleteDB()

    logging.info("Drawing graph in {}".format(args.out))

    if args.GaussianSmoothing:

        Graphics.plotDistributionWithGeneHistogram([x for x in range(-winBefore,winAfter+1-lenFilter)],lAveragePhases[0:(winBefore+winAfter+1-lenFilter)],lPhasesNb[0:(winBefore+winAfter+1-lenFilter)],lOtherGenesNb[0:(winBefore+winAfter+1-lenFilter)],out=args.out, title=args.title, xax=args.xax, yax=args.yax, yax2=args.zax)

    else:
        Graphics.plotDistributionWithGeneHistogram([x for x in range(-winBefore,winAfter+1)],lAveragePhases[0:(winBefore+winAfter+1)],lPhasesNb[0:(winBefore+winAfter+1)],lOtherGenesNb[0:(winBefore+winAfter+1)],out=args.out, title=args.title, xax=args.xax, yax=args.yax, yax2=args.zax)

