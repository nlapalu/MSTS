#!/usr/bin/env python

import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib import gridspec

class Graphics(object):

    def __init__(self):
        pass

    @staticmethod
    def plotDistribution(lXs, lYs, out="", title="", xax="", yax=""):
        """Draw a simple Distribution"""

        fig = plt.Figure(figsize=(20,20))
        fig.suptitle(title, fontsize=32)
        ax = fig.add_subplot(111)
        ax.plot(lXs,lYs)
        axis_font = {'size':'28'}
        ax.set_xlabel(xax, **axis_font)
        ax.set_ylabel(yax, **axis_font)
        ax.tick_params(labelsize=20)
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(out, dpi=80)

    @staticmethod
    def plotDistributionWithRegression(lXs, lYs, rXs, rYs, slope, intercept, R2, pval, meanPhase, stdvPhase,out="", title="", xax="", yax=""):
        """Draw a simple Distribution and add a regression subgraph"""

        fig = plt.Figure(figsize=(20,20))
        fig.suptitle(title, fontsize=32)
        ax = fig.add_subplot(111)
        ax.plot(lXs,lYs)
        annotation_font = {'size':16}
        for i in rYs:
            ax.annotate(str(i),xy=(i-10,lYs[i]+10), **annotation_font)
        ax.scatter([i for i in rYs],[lYs[i] for i in rYs],color='red',marker='+', s=200)
        axis_font = {'size':'28'}
        ax.set_xlabel(xax, **axis_font)
        ax.set_ylabel(yax, **axis_font)
        ax.tick_params(labelsize=20)
        ax.set_xlim(0,lXs[-1]+10)
        ax.set_ylim(min(lYs)-10,max(lYs)+10)
        # add regression graph
        ax2 = fig.add_axes([0.55,0.55,0.3,0.3])
       # ax2.axis([0,rXs[-1],0,rYs[-1]])
        ax2.set_xlim(0,rXs[-1]+2)
        ax2.scatter([x+1 for x in rXs],rYs)
        ax2.plot([x+1 for x in rXs], [x*slope+intercept for x in rXs], color='red')
        axis_font2 = {'size':'18'}
        ax2.set_title('linear regression', **axis_font2)
        ax2.set_xlabel('#peak', **axis_font2)
        ax2.set_ylabel('positions bp', **axis_font2)
        ax2.tick_params(labelsize=16)
        ax2.text(0.5,rYs[-1]*0.95 , 'Phase={}+-{}'.format(meanPhase,stdvPhase),
        verticalalignment='top', horizontalalignment='left',
        fontsize=15)
        ax2.text(0.5,rYs[-1]*0.9, 'R2={}'.format(R2),
        verticalalignment='top', horizontalalignment='left',
        fontsize=15)
        ax2.text(0.5,rYs[-1]*0.85, 'p-value={}'.format(pval),
        verticalalignment='top', horizontalalignment='left',
        fontsize=15)
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(out, dpi=80)


