#/usr/bin/env python

from __future__ import division
import sys
import numpy
import pyBigWig
import math
import argparse
import logging

from MSTS.version import __version__
from MSTS.Parser.ExpNucFileParser import ExpNucFileParser


from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fcluster
from scipy.stats import variation
from scipy import stats
import numpy as np

from scipy.stats import poisson
from MSTS.Graphics import Graphics

from statsmodels.sandbox.stats.multicomp import multipletests

import pandas as pd
import seaborn as sns; sns.set(color_codes=True)


#https://stackoverflow.com/questions/33944914/implementation-of-e-test-for-poisson-in-python/34013250#34013250

#https://stat.ethz.ch/R-manual/R-devel/library/stats/html/poisson.test.html

#https://stackoverflow.com/questions/33944914/implementation-of-e-test-for-poisson-in-python


# copied from statsmodels.stats.weightstats
def _zstat_generic2(value, std_diff, alternative):
    '''generic (normal) z-test to save typing

    can be used as ztest based on summary statistics
    '''
    zstat = value / std_diff
    if alternative in ['two-sided', '2-sided', '2s']:
        pvalue = stats.norm.sf(np.abs(zstat))*2
    elif alternative in ['larger', 'l']:
        pvalue = stats.norm.sf(zstat)
    elif alternative in ['smaller', 's']:
        pvalue = stats.norm.cdf(zstat)
    else:
        raise ValueError('invalid alternative')
    return zstat, pvalue


def poisson_twosample(count1, exposure1, count2, exposure2, ratio_null=1,
                      method='score', alternative='2-sided'):

    # shortcut names
    y1, n1, y2, n2 = count1, exposure1, count2, exposure2
    d = n2 / n1
    r = ratio_null
    r_d = r / d

    if method in ['score']:
        stat = (y1 - y2 * r_d) / np.sqrt((y1 + y2) * r_d)
        dist = 'normal'

    if dist == 'normal':
        return _zstat_generic2(stat, 1, alternative)


def quantileNormalization(lData):

    lSortedIndex = []
    lSortedValues = []
    for data in lData:
        sortedIndex = sorted(range(len(data)), key=lambda k: data[k])
        lSortedIndex.append(sortedIndex)
        sortedValues = [None]*len(data)
        for i,idx in enumerate(sortedIndex):
            sortedValues[i] = data[idx]
        lSortedValues.append(sortedValues)

    lSortedMeanValues = []
    for i,row in enumerate(lData[0]):
        lSortedMeanValues.append(numpy.mean([x[i] for x in lSortedValues]))

    lNormedValues = []
    for index in lSortedIndex:
        normedValues = [None]*len(index)
        for i,idx in enumerate(index):
            normedValues[idx] = lSortedMeanValues[i]
        lNormedValues.append(normedValues)

    return lNormedValues


def scallingNormalization(lData):

    lNormedValues = []
    lSum = []
    for i,data in enumerate(lData):
        lSum.append(0)
        for val in data:
            lSum[i] += val
    globalMean = np.mean(lSum)
    for i,data in enumerate(lData):
        lValues = []
        for val in data:
            lValues.append(val * globalMean/lSum[i])
        lNormedValues.append(lValues)

    return lNormedValues


if __name__ == "__main__":

    program = sys.argv[0]
    version = __version__
    description = 'Perform differential positioning analysis, \
                   under test'
    parser = argparse.ArgumentParser(prog=program)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--version', action='version', version='{} {}'.format(program,version))
    parser.add_argument("lNucFiles", help="File of files with labels, conditions, nucleosomes and bigwig files")
    parser.add_argument("-v", "--verbosity", type=int, choices=[1,2,3],
                        help="increase output verbosity 1=error, 2=info, 3=debug")
    parser.add_argument("-p", "--prefix", help="prefix for output files, default=[out.]", type=str, default="out")
    parser.add_argument("-a", "--alpha", help="alpha risk threshold to output results, default=0.05", type=float, default=0.05)
    parser.add_argument("-f", "--fc", help="Fold Change to output results, default=2.0", type=float, default=2.0)
    parser.add_argument("-t", "--threshold", help="Miminal sum of counts per condition to analyze position, default=5", type=int, default=5)
    parser.add_argument("-n", "--norm", help="Normalization method: quantile or scale, default=quantile", type=str, default="quantile")
    parser.add_argument("-y", "--nuctype", help="Nucleosome type to analyze, inclusive mode: very-well, well, fuzzy, bad, default=fuzzy", type=str, default="fuzzy")
    parser.add_argument("--bed", help="export differential positions in bed file", action="store_true", default=False)
    parser.add_argument("--merge", help="In case of several replicats with different peaks, positions could be merge in bed export to propose a differential area instead of several unique positions, default=False", action="store_true", default=False)

    args = parser.parse_args()

    logLevel='ERROR'
    if args.verbosity == 1:
        logLevel = 'ERROR'
    if args.verbosity == 2:
        logLevel = 'INFO'
    if args.verbosity == 3:
        logLevel = 'DEBUG'
    logger = logging.getLogger().setLevel(logLevel)

    if args.norm not in ["quantile", "scale"]:
        logging.error('Please choose the normalization method in this list: [quantile, scale]')
        sys.exit(0)

    if args.nuctype not in ["very-well", "well","fuzzy", "bad"]:
        logging.error('Please choose the nucleosome type to analyze in this list: [very-well, well, fuzzy, bad]')
        sys.exit(0)

    lNucTypes = ["very-well", "well","fuzzy", "bad"][0:["very-well", "well","fuzzy", "bad"].index(args.nuctype)+1]
    logging.info("Analyzing nucleosomes: {}".format(lNucTypes))

    prs = ExpNucFileParser(args.lNucFiles) 
    prs.parse()
    lFiles = prs.getlNucFiles()
    lWigFiles= prs.getlWigFiles()
    labels = prs.getlLabels()
    lConditions = prs.getlConditions()
    lPositions = set()


    # read all files
    for replicate in lFiles:
        #print replicate
        lnuc = []
        with open(replicate, 'r') as r:
            header = r.readline()
            for line in r:
                nuc = line.rstrip().split("\t")
                # extract intersting positions (center of defined conserved nuc)
                if nuc[7] in lNucTypes:
                    lPositions.add((nuc[0],int(nuc[1])+73))

    lData = []

#    ltmp = []
#    for i,val in enumerate(lPositions):
#        ltmp.append(val)
#        if i == 10000:
#            break
#    lPositions = ltmp 



    for replicate in lWigFiles:
        #print replicate
        bw = pyBigWig.open(replicate)
        values = []
        for position in lPositions:
            values.append(bw.values(position[0], position[1]-1, position[1])[0])
        lData.append(values)

    lNormData = []
    if args.norm == "quantile":
        lNormData = quantileNormalization(lData)
    elif args.norm == "scale":
        lNormData = scallingNormalization(lData)
    else:
        raise  Exception("Problem unknown normalization method")

    logging.info("hierarchical clustering of raw data: {}".format("{}_diff_hc_raw.png".format(args.prefix)))
    x = np.array(lData)
    Z = linkage(x, 'ward')
    Graphics.plotHierarchicalClustering(Z, 0.3*max(Z[:,2]), labels=labels, out="{}_diff_hc_raw.png".format(args.prefix), title="Hierarchical clustering of raw data for {}".format(args.prefix), xax="replicats", yax="")

    logging.info("hierarchical clustering of normalized data: {}".format("{}_diff_hc_norm.png".format(args.prefix)))
    x = np.array(lNormData)
    Z = linkage(x, 'ward')
    Graphics.plotHierarchicalClustering(Z, 0.3*max(Z[:,2]), labels=labels, out="{}_diff_hc_norm.png".format(args.prefix), title="Hierarchical clustering of normalized data for {}".format(args.prefix), xax="replicats", yax="")

    n1 = 0
    n2 = 0
    for i,pos in enumerate(lPositions):
        analyze = True
        for condition in lConditions:
            if sum([lData[rep][i] for rep in condition]) < args.threshold:
                analyze = False
        if analyze:
            n1 += sum([lNormData[rep][i] for rep in lConditions[0]])
            n2 += sum([lNormData[rep][i] for rep in lConditions[1]])

    lFinal = []
    lpvalues = []

    for i,pos in enumerate(lPositions):
        analyze = True
        for condition in lConditions:
            if sum([lData[rep][i] for rep in condition]) < args.threshold:
                analyze = False

        if analyze:
            c1 = sum([lNormData[rep][i] for rep in lConditions[0]])
            c2 = sum([lNormData[rep][i] for rep in lConditions[1]])
            nb1 = len(lConditions[0])
            nb2 = len(lConditions[1])
            mean1 = c1/nb1
            mean2 = c2/nb2
            var1 = variation([lNormData[rep][i] for rep in lConditions[0]]) * 100
            var2 = variation([lNormData[rep][i] for rep in lConditions[1]]) * 100
            FC = mean2/mean1
            s, pv = poisson_twosample(c1,n1,c2,n2,method='score')
            #pv = poisson.pmf(round(min(lMeanNormData[0][i],lMeanNormData[1][i])),round(max(lMeanNormData[0][i],lMeanNormData[1][i])))
    
            lFinal.append([pos[0],pos[1],mean1,var1,mean2,var2,FC,pv])
            lpvalues.append(pv)
    lFDRvalues = multipletests(lpvalues, alpha=args.alpha, method='fdr_bh')[1]
    for i,entry in enumerate(lFinal):
        entry.append(lFDRvalues[i])
    FCSortedlFinal = sorted(lFinal, key=lambda k: (k[6]))

    logging.info("Histogram of p-values in {}".format("{}_hist_pvalues.png".format(args.prefix)))
    Graphics.plotPValHistogram(lpvalues,[ x*0.05 for x in range(0,21)], out="{}_hist_pvalues.png".format(args.prefix), title="Histogram of p-values for {}".format(args.prefix), xax="p-values", yax="positions")


    FCRetainedlFinal = []
    for entry in FCSortedlFinal:
        if entry[8] < args.alpha and (entry[6] < 1/args.fc or entry[6] > args.fc):
            FCRetainedlFinal.append(entry)

    logging.info("Results exported in {}.diff".format(args.prefix))
    with open("{}.diff".format(args.prefix), "w") as f:
        f.write("#ref\tpos\tmean1\tCV1\tmean2\tCV2\tFC\tpval\tpaj\n")
        for entry in FCRetainedlFinal:
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(entry[0],entry[1],"{0:1.1f}".format(entry[2]),"{0:1.1f}".format(entry[3]),"{0:1.1f}".format(entry[4]), "{0:1.1f}".format(entry[5]),"{0:1.2f}".format(entry[6]),"{0:1.2e}".format(entry[7]), "{0:1.2e}".format(entry[8])))
    f.close()

    logging.info("{} positions tested".format(len(FCSortedlFinal)))
    logging.info("{} positions detected with a potential difference in signal".format(len(FCRetainedlFinal)))

    if args.bed:
        logging.info("Results in bed format exported in {}.bed".format(args.prefix))
        with open("{}.bed".format(args.prefix), 'w') as f:
            if args.merge:
                PosSortedlFinal = sorted(FCRetainedlFinal, key=lambda k: (k[0],k[1]))
                currentPos = (PosSortedlFinal[0][0],PosSortedlFinal[0][1]-1,PosSortedlFinal[0][1])
                for pos in PosSortedlFinal[1:]:
                    if pos[0] == currentPos[0] and (pos[1]-1 >= currentPos[1]-30 and pos[1] <= currentPos[2]+30):
                        logging.info("merge: {} with {}".format(currentPos,(pos[0],pos[1])))
                        currentPos = (currentPos[0],min(currentPos[1], pos[1]-1),max(currentPos[2], pos[1]))
                    else:
                        f.write("{}\t{}\t{}\n".format(currentPos[0],currentPos[1],currentPos[2]))
                        currentPos = (pos[0],pos[1]-1,pos[1])
            else:
                for entry in sorted(FCRetainedlFinal, key=lambda k: (k[0],k[1])):
                    f.write("{}\t{}\t{}\n".format(entry[0],entry[1]-1,entry[1]))
        f.close()
