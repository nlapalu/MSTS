#!/usr/bin/env python

import sys
import logging
import time
import argparse
import pyBigWig
import math
import random
import copy

from collections import Counter

from fastcluster import *

from matplotlib import pyplot as plt
#from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fcluster
from scipy.stats import variation
from scipy import stats
import numpy as np

from sklearn.cluster import KMeans

from MSTS.version import __version__
from MSTS.Parser.SimpleGffParser import SimpleGffParser
from MSTS.Db.FeatureDB import FeatureDB
from MSTS.Graphics import Graphics


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
        ltmp = [0.0]*(windowWidth*stdev)
        data.extend(ltmp)
        ltmp.extend(data)
        data = ltmp

        for i in range(0,len(smoothed)):
            for j in range(0,len(filter)):
                smoothed[i] += data[i + j] * filter[j]
        smoothed[0:windowWidth*stdev] = [None]*(windowWidth*stdev)
        smoothed[-(windowWidth*stdev):] = [None]*(windowWidth*stdev)

        return smoothed
 

def prepare_export(fOut):
    """Prepare export, remove file and add header"""

    with open(fOut, 'w') as f:
        f.write("seq\tstart\tend\tmean\tstdev\tpeak\tcluster\tpositioning\n")
    f.close()


def export(fOut,lNucleosomes):
    """Write nucleosomes in out file"""

    with open(fOut,'a') as f:
        for nuc in lNucleosomes:
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(nuc[0],nuc[1]+1,nuc[2],nuc[3],nuc[4],nuc[5],nuc[7],nuc[8]))
    f.close()

def exportBed(fOut,lNucleosomes):
    """Write nucleosomes in bed file"""

    with open(fOut,'w') as f:
        for nuc in lNucleosomes:
            f.write("{}\t{}\t{}\n".format(nuc[0],nuc[1],nuc[2]))
    f.close()
 
  
def printWigHeader():
    """write wig header"""

    try:
        with open(args.prefix + ".wig", 'w') as f:
            f.write("track type=wiggle_0 name=fileName description=fileName\n")
        f.close()
    except:
        logging.error("Can not write results in {}.wig".format(args.prefix))

def printWig(seq, lcov):
    """write output in wig format"""

    try:
        with open(args.prefix + ".wig", 'a') as f:
            f.write("variableStep chrom={}\n".format(seq))
#            if self.keepPosBedFile:
#                for pos,cov in enumerate(lcov):
#                    if self.isFragmentFullyIncludeInABedTrack(seq,pos,pos+1):
#                        f.write("{}\t{}\n".format(pos+1,cov))
#            else:
            for pos,cov in enumerate(lcov):
                if cov != None:
                    f.write("{}\t{}\n".format(pos+1,cov))
        f.close()
    except:
        logging.error("Can not write results in {}.wig".format(args.prefix))


def clusterClassifier(lStats):
    """Classify clusters"""

    lKClassif = ['NA']*len(lStats)
    for i,k in enumerate(lStats):
       # pseudo-area todo with AUC
        tot = sum([k[0][j] for j in range(26,121)])
        if tot >= 0.85:
            lKClassif[i] = 'very-well'
        elif tot >= 0.80:
            lKClassif[i] = 'well'
        elif tot >= 0.65:
            lKClassif[i] = 'fuzzy'
        else:
            lKClassif[i] = 'bad'
        # middle position > 60 and < 86
        peak = 0.0
        peakPos = 0
        for j,mean in enumerate(k[0]):
            #print mean
            if mean > peak:
                peakPos = j
                peak = mean
        if peakPos < 60 or peakPos > 86:
            lKClassif[i] = 'bad'
            
    return lKClassif

def clusterClassifierRelaxed(lStats):
    """Classify clusters"""

    lKClassif = ['NA']*len(lStats)
    for i,k in enumerate(lStats):
       # pseudo-area todo with AUC
        tot = sum([k[0][j] for j in range(26,121)])
        if tot >= 0.85:
            lKClassif[i] = 'very-well'
        elif tot >= 0.75:
            lKClassif[i] = 'well'
        elif tot >= 0.65:
            lKClassif[i] = 'fuzzy'
        else:
            lKClassif[i] = 'bad'
        # middle position > 60 and < 86
        peak = 0.0
        peakPos = 0
        for j,mean in enumerate(k[0]):
            #print mean
            if mean > peak:
                peakPos = j
                peak = mean
        if peakPos < 60 or peakPos > 86:
            lKClassif[i] = 'bad'
            
    return lKClassif
        

if __name__ == "__main__":

    program = sys.argv[0]
    version = __version__
    description = 'This tool detects and classifies nucleosomes from a bigWig \
                   file. First, It performs a gaussian smoothing on the data. Then It launches two \
                   clustering steps: 1/ a hierarchical clustering, done in several iterations with \
                   random subsets of data. This step defined the most probable number of clusters \
                   expected 2/ a KMeans clustering on all data with the previously defined number \
                   of cluster. \
                   Clusters are then exported in several formats (bed, wig).'

    parser = argparse.ArgumentParser(prog=program)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--version', action='version', version='{} {}'.format(program,version))

    parser.add_argument("bigWig", help="Input bigWig File", type=str)

    parser.add_argument("-w","--windowWidth", help="window size for Gaussian smoothing, default=3", type=int, default=3)
    parser.add_argument("-sd","--stdev", help="stdev for Gaussian smoothing, default=20", type=int, default=20)
    parser.add_argument("-b", "--bigBed", help="bigBedFile, use to limit detections to specific regions", type=str, default=None)
    parser.add_argument("-p", "--prefix", help="prefix for output files, default=[out.]",type=str, default="out")
    parser.add_argument("--wig", help="output wig file", action="store_true", default=False)
    parser.add_argument("--bed", help="output bed files per cluster", action="store_true", default=False)
    parser.add_argument("-df","--distanceFactor", help="factor used to compute distance when iterating with hierarchical clustering, (factor_distance X max(distances)  default=0.3, range[0.1-0.9]", type=float, default=0.3)
    parser.add_argument("-nbi","--nbIterations", help="nb iteration for hierarchical clustering, default=500", type=int, default=500)
    parser.add_argument("-nbn","--nbNucHC", help="nb nucleosomes to sample for hierarchical clustering, default=2000", type=int, default=2000)
    parser.add_argument("-mic", "--minCov", help="minimum coverage of dyad to keep the nucleosome, by default 20% of median, default=0.2", type=float, default=0.2) 
    parser.add_argument("-ov", "--overlap", help="Allow overlap between nucleosomes, default=30", type=int, default=30)
    parser.add_argument("--refine", help="Refine detection on nucleosome classified as fuzzy and bad, default=true", action="store_true", default=False)
    parser.add_argument("-t", "--title", help="title text", type=str, default="Nucleosome profil clustering")
    parser.add_argument("-x", "--xax", help="x axis text", type=str, default="window bp")
    parser.add_argument("-y", "--yax", help="y axis text", type=str, default="proportion, %")
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

    if args.distanceFactor < 0.1 or args.distanceFactor > 0.9:
        logging.error('DistanceFactor out of authorized range: [0.1-0.9]')


    prepare_export("{}.nucleosomes.txt".format(args.prefix))
    if args.refine:
        prepare_export("{}.refine.nucleosomes.txt".format(args.prefix))

    if args.wig:
        printWigHeader()

    buffSize = 20000000
#    buffSize = 200000
    logging.info('Buffer size: {} bases'.format(buffSize))

    bw = pyBigWig.open(args.bigWig)

    lAllProfils = []
    lAllProfilNorms = []
    lAllNucleosomes = []
    for chrom in bw.chroms():
        lNucleosomes = []
        lProfils = []
        lProfilNorms = []
        lChrSmoothedValues = []
        logging.info('Reading sequence: {}'.format(chrom))
        if args.bigBed:
            bb = pyBigWig.open(args.bigBed)
            logging.error("Not yet implemented")
            sys.exit(1)
            if chrom not in bb.chroms().keys():
                continue
                #TODO bigbed 
        else:
            start = 0
            nb_chunks = 0
            while (start < bw.chroms(chrom)): 
                stop = min(start+buffSize, bw.chroms(chrom))
                logging.info('Requesting values: {}:{}-{}'.format(chrom,start,stop))
                values = bw.values(chrom, start, stop)
                #print len(values)
                #print values[0]

                #TODO Handle nan values
                #for idx,val in enumerate(values):
                #    if math.isnan(val):
         

                logging.info("Smoothing data with Gaussian blur, window: {}, stdev: {}".format(args.windowWidth, args.stdev))
                lSmoothedValues = gaussianSmoothing(values, args.windowWidth, args.stdev)
                lChrSmoothedValues.extend(lSmoothedValues)
                
                # Fix python 3 ; sort error with NoneType
                #lSortedIndices = sorted(range(len(lSmoothedValues)), key=lSmoothedValues.__getitem__, reverse=True)
                lModifiedValues = []
                for v in lSmoothedValues :
                    if v != None:
                        lModifiedValues.append(v)
                    else:
                        lModifiedValues.append(-999)
                lSortedIndices = sorted(range(len(lModifiedValues)), key=lModifiedValues.__getitem__, reverse=True)

                for idx in lSortedIndices:
                    
                   # Fix python 3 ; compare error with NoneType
                   # if lSmoothedValues[idx] > 0 and values[idx] != 0:
                    if lModifiedValues[idx] > 0 and values[idx] != 0:
                        #print max(idx-73+start,start)
                        #print min(idx+73+1+start,stop)
                        mean = np.mean([values[i] for i in range(max(idx-73+start,start)-buffSize*nb_chunks,min(idx+73+1+start,stop)-buffSize*nb_chunks)])
                        stdev = np.std([values[i] for i in range(max(idx-73+start,start)-buffSize*nb_chunks,min(idx+73+1+start,stop)-buffSize*nb_chunks)])
                        nuc = [chrom,max(idx-73+start,start),min(idx+73+1+start,stop), mean, stdev,values[idx],[]]
                        lNucleosomes.append(nuc)
                        lAllNucleosomes.append(nuc)
                        lProfils.append([values[x] for x in range(idx-73,idx+73+1)])
                        lAllProfils.append([values[x] for x in range(idx-73,idx+73+1)])
                        tot = sum([values[x] for x in range(idx-73,idx+73+1+1)])
                        lProfilNorms.append([values[x]/tot for x in range(idx-73,idx+73+1)])
                        lAllProfilNorms.append([values[x]/tot for x in range(idx-73,idx+73+1)])
                        nuc[6] = [values[x]/tot for x in range(idx-73,idx+73+1)]
                        limit1 = max(idx-147+args.overlap,0)
                        limit2 = min(idx+147-args.overlap,len(lSmoothedValues))
                        for i in range(limit1,limit2):
                            lSmoothedValues[i] = 0 

                start = stop + 1
                nb_chunks += 1

            logging.info("{} nucleosomes detected and positioned for {}".format(len(lNucleosomes),chrom))

            if args.wig:
                printWig(chrom,lChrSmoothedValues)

##            lSortedNucleosomes = sorted(lNucleosomes, key=lambda nuc: nuc[1])
##            export("{}.nucleosomes.txt".format(args.prefix),lSortedNucleosomes)
    logging.info("{} nucleosomes detected ans positioned for all sequences".format(len(lAllProfilNorms)))

    median = np.median([nuc[5] for nuc in lAllNucleosomes])

    lCleanNuc = []
    lCleanProfils = []
    lCleanNormProfils = []
    nbNucFiltered = 0  
    for i,nuc in enumerate(lAllNucleosomes):
        if nuc[5] > median*args.minCov:
            lCleanNuc.append(nuc)
            lCleanProfils.append(lAllProfils[i])
            lCleanNormProfils.append(lAllProfilNorms[i])
        else:
            nbNucFiltered += 1

    lAllNucleosomes = lCleanNuc
    lAllProfils = lCleanProfils
    lAllProfilNorms = lCleanNormProfils
            
    logging.info("{} nucleosomes detected ans positioned for all sequences after filter of mincov: {}, median: {}, {} nucleosomes removed".format(len(lAllProfilNorms),median*args.minCov,median,nbNucFiltered))

 

    # step 1: hierarchical clustering on subset (X iter)
    lhClusters = []
    logging.info("hierarchical clustering: {} iterations with {} nucleosomes".format(args.nbIterations, args.nbNucHC))
    for nbIt in range(0,args.nbIterations):
        lHProfils = []
        for i in range(0,args.nbNucHC):
            lHProfils.append(lAllProfilNorms[random.randint(0,len(lAllProfilNorms)-1)])
        x = np.array(lHProfils)
        Z = linkage(x, 'ward')
        t = args.distanceFactor*max(Z[:,2])
        ar = fcluster(Z,t,criterion='distance')

        hClusters = set()
        for i in np.nditer(ar):
           hClusters.add(int(i))
        lhClusters.append(len(hClusters))
        logging.debug("hierachical clustering iter: {}, nb clusters: {}".format(nbIt,len(hClusters)))

    dClusters = Counter(lhClusters)
    nbClusters = sorted(dClusters.items(), key=lambda x:x[1], reverse=True)[0][0]
    logging.info("{} clusters defined after {} iterations".format(nbClusters,args.nbIterations))
    if nbClusters > 25:
        nbClusters = 25
        logging.info("Too many clusters, number limited to {} clusters".format(nbClusters))

    # step 2: k-means clustering with pre-defined nb clusters
    logging.info("Performing K-Means clustering") 
    lkmeans = KMeans(n_clusters=nbClusters,random_state=185).fit_predict(lAllProfilNorms)
    lStats = []
    lKs = [0]*nbClusters
    lAllkNucleosomes = []
    logging.info("Extracting info from clustering")
    for k in range(0,nbClusters):
        lkProfils = []
        lkNucleosomes = []
        j = 0
        for i in np.nditer(lkmeans):
            if i == k:
                lkProfils.append(lAllProfilNorms[j])
                lkNucleosomes.append(lAllNucleosomes[j])
                # add cluster info
                lAllNucleosomes[j].append(k)
                lAllkNucleosomes.append(lAllNucleosomes[j])
                lKs[k] += 1
            j+=1 

        lmeank = [0.0]*147
        lmink = [0.0]*147
        lmaxk = [0.0]*147
        lquartile1 = [0.0]*147
        lquartile3 = [0.0]*147
        for i in range(0,147):
            lmeank[i] = np.mean([x[i] for x in lkProfils])
            lmink[i] = min([x[i] for x in lkProfils])
            lmaxk[i] = max([x[i] for x in lkProfils])
            lquartile1[i] = np.percentile([x[i] for x in lkProfils],25)
            lquartile3[i] = np.percentile([x[i] for x in lkProfils],75)
        lStats.append([lmeank,lmink,lmaxk, lquartile1, lquartile3])


    lKClassif = clusterClassifier(lStats)
    # add position info
    for nuc in lAllkNucleosomes:
        nuc.append(lKClassif[nuc[7]])

    lSortedAllkNucleosomes = sorted(lAllkNucleosomes, key=lambda nuc: (nuc[0],nuc[1]))
#    export("k{}.clusters.txt".format(idx),lSortedlkNucleosomes)
    if not args.refine:
        export("{}.nucleosomes.txt".format(args.prefix),lSortedAllkNucleosomes)
    export("{}.nucleosomes.txt".format(args.prefix),lSortedAllkNucleosomes)
#        lSortedlkNucleosomes = sorted(lkNucleosomes, key=lambda nuc: (nuc[1],nuc[2]))
#        export("k{}.clusters.txt".format(idx),lSortedlkNucleosomes)
#    lClassification = clusterClassifier(lKs)

    if args.bed:
        if not args.refine:
            for k in range(0,nbClusters):
                exportBed("{}.k{}-{}.cluster.bed".format(args.prefix,k,lKClassif[k]),[x for x in lSortedAllkNucleosomes if x[7] == k])

    logging.info("Drawing clusters in {}.clusters.png".format(args.prefix))
    Graphics.plotDistributionWithLimits([i for i in range(1,148)],lStats,lKClassif,out="{}.clusters.png".format(args.prefix), title=args.title,xax=args.xax, yax=args.yax,legend=lKs)


    if args.refine:
        lStatInit = []
        lKClassifInit = []
        lKInit = []
        dReorder = {}
        Idx = 0
        for i,val in enumerate(lStats):
            if lKClassif[i] not in ["fuzzy", "bad"]:
                lStatInit.append(lStats[i])
                lKClassifInit.append(lKClassif[i])
                lKInit.append(lKs[i])
                dReorder[i] = Idx
                Idx += 1

        lRefinedNuc = copy.deepcopy([ nuc for nuc in lSortedAllkNucleosomes if nuc[8] not in ["fuzzy", "bad"] ])
        for nuc in lRefinedNuc:
            nuc[7] = dReorder[nuc[7]] 
        #print len(lRefinedNuc)

        NbClusInit = len(lKInit)

        for Kclass in ["fuzzy", "bad"]:
            
            logging.info("refine {} nucleosomes: re-clustering".format(Kclass))

            lAllNucleosomesFuzzy = copy.deepcopy([nuc for nuc in lSortedAllkNucleosomes if nuc[8] == Kclass ])
            if len(lAllNucleosomesFuzzy) == 0:
                continue
            lnucFuzzy = [nuc[6] for nuc in lSortedAllkNucleosomes if nuc[8] == Kclass ] 
        #    print len(lAllNucleosomesFuzzy) 

            lhClusters = []
            logging.info("hierarchical clustering: {} iterations with {} nucleosomes".format(args.nbIterations, args.nbNucHC))
            for nbIt in range(0,args.nbIterations):
                lHProfils = []
                for i in range(0,args.nbNucHC):
                    lHProfils.append(lnucFuzzy[random.randint(0,len(lnucFuzzy)-1)])
                x = np.array(lHProfils)
                Z = linkage(x, 'ward')
                t = args.distanceFactor*max(Z[:,2])
                ar = fcluster(Z,t,criterion='distance')

                hClusters = set()
                for i in np.nditer(ar):
                   hClusters.add(int(i))
                lhClusters.append(len(hClusters))
                logging.debug("hierachical clustering iter: {}, nb clusters: {}".format(nbIt,len(hClusters)))

            dClusters = Counter(lhClusters)
            nbClusters = sorted(dClusters.items(), key=lambda x:x[1], reverse=True)[0][0]
            logging.info("{} clusters defined after {} iterations".format(nbClusters,args.nbIterations))
            if nbClusters > 25:
                nbClusters = 25
                logging.info("Too many clusters, number limited to {} clusters".format(nbClusters))



            # step 2: k-means clustering with pre-defined nb clusters
            logging.info("Performing K-Means clustering") 
            lkmeans = KMeans(n_clusters=nbClusters,random_state=185).fit_predict(lnucFuzzy)
            lStats = []
            lKs = [0]*nbClusters
            lAllkNucleosomesFuzzy = []
            logging.info("Extracting info from clustering")
            for k in range(0,nbClusters):
                lkProfils = []
               # lkNucleosomes = []
                j = 0
                for i in np.nditer(lkmeans):
                    if i == k:
                        lkProfils.append(lnucFuzzy[j])
                #        lkNucleosomes.append(lAllNucleosomesFuzzy[j])
                        # add cluster info
                        lAllNucleosomesFuzzy[j][7] = k + NbClusInit
                        lAllkNucleosomesFuzzy.append(lAllNucleosomesFuzzy[j])
                        #if lAllNucleosomesFuzzy[j][1] == 28174:
                        #    print "ye"
                        lKs[k] += 1
                    j+=1 

                lmeank = [0.0]*147
                lmink = [0.0]*147
                lmaxk = [0.0]*147
                lquartile1 = [0.0]*147
                lquartile3 = [0.0]*147
                for i in range(0,147):
                    lmeank[i] = np.mean([x[i] for x in lkProfils])
                    lmink[i] = min([x[i] for x in lkProfils])
                    lmaxk[i] = max([x[i] for x in lkProfils])
                    lquartile1[i] = np.percentile([x[i] for x in lkProfils],25)
                    lquartile3[i] = np.percentile([x[i] for x in lkProfils],75)
                lStats.append([lmeank,lmink,lmaxk, lquartile1, lquartile3])


            lKClassif = clusterClassifierRelaxed(lStats)
            # add position info
            for nuc in lAllkNucleosomesFuzzy:
                nuc[8] = lKClassif[nuc[7] - NbClusInit]




#        lSortedAllkNucleosomesFuzzy = sorted(lAllkNucleosomesFuzzy, key=lambda nuc: (nuc[0],nuc[1]))
#    export("k{}.clusters.txt".format(idx),lSortedlkNucleosomes)
#    export("{}.nucleosomes.txt".format(args.prefix),lSortedAllkNucleosomes)
#        lSortedlkNucleosomes = sorted(lkNucleosomes, key=lambda nuc: (nuc[1],nuc[2]))
#        export("k{}.clusters.txt".format(idx),lSortedlkNucleosomes)
#    lClassification = clusterClassifier(lKs)

            logging.info("Drawing clusters in {}.clusters.{}.png".format(args.prefix, Kclass))
            Graphics.plotDistributionWithLimits([i for i in range(1,148)],lStats,lKClassif,out="{}.clusters.refine.{}.png".format(args.prefix,Kclass), title="{} - refine {}".format(args.title,Kclass),xax=args.xax, yax=args.yax,legend=lKs)


            NbClusInit += len(lKs)
            lKInit.extend(lKs)
            lRefinedNuc.extend(lAllkNucleosomesFuzzy)
            lStatInit.extend(lStats)
            lKClassifInit.extend(lKClassif)
           

          
#        lSortedAllkNucleosomesNew = []
        lSortedRefinedNuc = sorted(lRefinedNuc, key=lambda nuc: (nuc[0],nuc[1]))
#        for i,nuc in enumerate(lSortedAllkNucleosomes):
#            for j,nuc2 in enumerate(lSortedRefinedNuc[i:]):
#                if (nuc[0],nuc[1]) == (nuc2[0],nuc2[1]):
#                    lSortedAllkNucleosomesNew.append(nuc2)
#                    break
#                else:
#                    lSortedAllkNucleosomesNew.append(nuc)
#                    break
#        print "BUG EXPORT"
        export("{}.refine.nucleosomes.txt".format(args.prefix),lSortedRefinedNuc)


        logging.info("Drawing refined clusters in {}.clusters.refine.png".format(args.prefix))
        Graphics.plotDistributionWithLimitsRefine([i for i in range(1,148)],lStatInit,lKClassifInit,out="{}.clusters.refine.png".format(args.prefix), title="{} - refined clusters".format(args.title),xax=args.xax, yax=args.yax,legend=lKInit)

        if args.bed:
            for k in range(0,NbClusInit):
                exportBed("{}.k{}-{}.cluster.refine.bed".format(args.prefix,k,lKClassifInit[k]),[x for x in lSortedRefinedNuc if x[7] == k])


