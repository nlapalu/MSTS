#!/usr/bin/env python

import logging
import sys
import argparse

from MSTS.version import __version__
from MSTS.TPMCounter import TPMCounter

if __name__ == '__main__':

    program = sys.argv[0]
    version = __version__
    description = 'Count TPM for each transcript from a BAM file It requires \
                   an annotation file in gff format, a sorted/indexed BAM file\
                   The program computes the mean fragment length for each transcript\
                   and the mean for all. Then, it computes the effective\
                   length and the TPM'
    parser = argparse.ArgumentParser(prog=program)
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("--version", action='version', version='{} {}'.format(program,version))
    parser.add_argument("BamFile", help="BamFile", type=str)
    parser.add_argument("AnnotFile", help="Annotation file in gff3 format", type=str)

    parser.add_argument("-t","--featType", help="Feature type choice for counts [exon,cds], default=exon", type=str, default="exon")
    parser.add_argument("-m","--minNbFrags", help="Minimum number of fragment per trancript to \
                        compute mean fragment length, [default=3]", type=int,
                        default=3)
    parser.add_argument("-s","--stranded", help="Define which fragments take into account\
                        for computing mean fragment length and counts. If your reads are \
                        not oriented, set to \"no\" as default, if not use \"yes\" or \
                        \"reverse\". [default=no]", type=str, default="no")
    #parser.add_argument("-c","--countFile", help="file with counts (as Htseq-count, ...) \
    #                    instead of internal count computation. Format: at least 2 columns \
    #                    with no header")
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


    if args.featType.upper() not in ["EXON","CDS"]:
        logging.error("available featType: exon or cds, please change your parameter")
        sys.exit(1)

    i = TPMCounter(args.BamFile, args.AnnotFile, logLevel)
    #i.run(args.minNbFrags, args.stranded, countFile=args.countFile)
    i.run(args.minNbFrags, args.stranded, countFile=None, featType=args.featType.upper())
