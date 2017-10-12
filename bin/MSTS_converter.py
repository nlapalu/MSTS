#/usr/bin/env python

import sys
import logging
import argparse

from MSTS.version import __version__
from MSTS.BamConverter import BamConverter

if __name__ == "__main__":

    program = sys.argv[0]
    version = __version__
    description = 'Convert Bam file to wig/bed/cov/size file, \
                   Your Bam file must be sorted and indexed, \
                   You can use samtools sort / samtools index'

    parser = argparse.ArgumentParser(prog=program)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--version', action='version', version='{} {}'.format(program,version))
    parser.add_argument("BamFile", help="Input Bam file", type=str)
    parser.add_argument("-p", "--prefix", help="prefix for output files, default=[out.]",
                        type=str, default="out")
    parser.add_argument("-m","--mode", help="counting mode type: [single, fragment, fragment-middle], default=single", type=str, default="single")
    parser.add_argument("--bed", help="output bed file", action="store_true", default=False)
    parser.add_argument("--wig", help="output wig file", action="store_true", default=False)
    parser.add_argument("--cov", help="output cov file", action="store_true", default=False)
    parser.add_argument("--size", help="output size file", action="store_true", default=False)
    parser.add_argument("-g", "--genome", help="Genome file used to sort references, structure: <refname><TAB><size>", type=str, default="")
    parser.add_argument("-w", "--window", help="For fragment middle mode +- window size, default=0", type=int, default=0)
    parser.add_argument("-k","--keepPosBedFile", help="bed file used to filter positions, keep position in file, warning only usable with seq < 100Mb", type=str)
    parser.add_argument("--minFSize", help="minimum fragment size to keep", type=int, default=0)
    parser.add_argument("--maxFSize", help="maximum fragment size to keep", type=int, default=2000)
    parser.add_argument("--minDepCov", help="minimum depth of coverage to keep, warning a posteriori filter, compute on reduce fragment if mode=fragment-middle", type=int, default=0)
    parser.add_argument("--maxDepCov", help="maximum depth of coverage to keep, warning a posteriori filter, compute on reduce fragment if mode=fragment-middle", type=int, default=1000000)
#    parser.add_argument("--primary-only", help="Only primary mappings are considered",
 #                       action="store_true")
   ## add only mapped reads !
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

    if args.mode not in ["single", "fragment", "fragment-middle"]:
        logging.error('Please choose the counting mode in this list: [single, fragment, fragment-middle]')
        sys.exit(0)
    else:
        if args.bed or args.wig or args.cov or args.size:
            i = BamConverter(args.BamFile, args.mode, args.prefix, args.bed, args.wig, args.cov, args.size, args.genome, args.window, args.keepPosBedFile, args.minFSize, args.maxFSize, args.minDepCov, args.maxDepCov, logLevel)
            i.convert()
        else:
            logging.error('Please choose at least one output format: [--bed, --wig, --cov, --size]')
            sys.exit(0)
