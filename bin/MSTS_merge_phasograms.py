#!/usr/bin/env python

import sys
import re
import argparse
import logging

from MSTS.version import __version__
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg


if __name__ == "__main__":

    program = sys.argv[0]
    version = __version__
    description = 'Merge phasograms from list'

    parser = argparse.ArgumentParser(prog=program)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--version', action='version', version='{} {}'.format(program,version))
    parser.add_argument("fileList", help="Input file list of phasograms to merge, format:<file>\tab<legend>", type=str)
    parser.add_argument('--colors', nargs='*', default=["red","blue","green","orange","purple","maroon"], help='Specify your colors instead of the 6 default values: ["red","blue","green","orange","purple","maroon"]', type=str)
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


    llXs = []
    llYs = []
    legend = []
    with open(args.fileList) as f:
        for line in f:
            if not re.match('^#', line):
                if line.rstrip() == '' :
                    continue
                values = line.rstrip().split('\t')
                if (len(values) != 2) :
                    logging.error("file of file: {} badly formatted".format(args.fileList))
                    sys.exit(1)
                legend.append(values[1])
                try:
                    with open(values[0]) as phaso:
                        lYs = []
                        lXs = []
                        for entry in phaso:
                            values2 = entry.rstrip().split('\t')
                            if values2[1] != "None":
                                lYs.append(float(values2[1]))
                            else:
                                lYs.append(None)
                            lXs.append(int(values2[0]))
                        llYs.append(lYs)
                    if llXs:
                        if llXs != lXs:
                            print("Error indice list")
                            sys.exit(1)
                    else:
                        llXs = lXs
                except:
                     logging.error("Error when reading {}, missing file or format error".format(values[0]))


    fig = plt.Figure(figsize=(20,20))
    fig.suptitle("merged phasogram", fontsize=32)
    ax = fig.add_subplot(111)
    for i,val in enumerate(llYs):
        ax.plot(lXs,val,color=args.colors[i%len(args.colors)])
    axis_font = {'size':'28'}
    ax.set_xlabel("window, bp", **axis_font)
    ax.set_ylabel("signal coverage", **axis_font)
    ax.tick_params(labelsize=20)
    ax.legend(legend, fontsize=20)
    canvas = FigureCanvasAgg(fig)
    canvas.print_figure("phaso.merge.png", dpi=80)
