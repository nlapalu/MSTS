#!/usr/bin/env python

import math
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib import gridspec

from scipy.cluster.hierarchy import dendrogram
from scipy.stats import variation

class Graphics(object):


    lColors = ['palevioletred','darkorchid','royalblue','darkturquoise','mediumspringgreen','olivedrab','gold','sandybrown','red','silver','m','blue','lightseagreen','chartreuse','darkkhaki','bisque','sienna','firebrick','gray','plum','darkcyan','darkgreen','orange','lightcoral','yellow']

    def __init__(self):
        pass

    @staticmethod
    def plotDistribution(lXs, lYs, out="", title="", xax="", yax="", color="blue", legend="", grid=[]):
        """Draw a simple Distribution"""

        fig = plt.Figure(figsize=(20,20))
        fig.suptitle(title, fontsize=32)
        ax = fig.add_subplot(111)
        ax.plot(lXs,lYs, color=color)
        if legend:
            ax.legend(legend, fontsize=22)
        for line in grid:
            ax.axvline(x=line, linestyle='dashed', linewidth=1, color='black')
        axis_font = {'size':'28'}
        ax.set_xlabel(xax, **axis_font)
        ax.set_ylabel(yax, **axis_font)
        ax.tick_params(labelsize=20)
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(out, dpi=80)


    @staticmethod
    def plotMultiDistribution(lXs, llYs, out="", title="", xax="", yax="", color=['blue','green'],legend="", grid=[]):
        """Draw multi Distributions"""

        fig = plt.Figure(figsize=(20,20))
        fig.suptitle(title, fontsize=32)
        ax = fig.add_subplot(111)
        for i,val in enumerate(llYs):
            ax.plot(lXs,val,color=color[i])
        axis_font = {'size':'28'}
        ax.set_xlabel(xax, **axis_font)
        ax.set_ylabel(yax, **axis_font)
        ax.tick_params(labelsize=20)
        if legend:
            ax.legend(legend, fontsize=22)
        for line in grid:
            ax.axvline(x=line, linestyle='dashed', linewidth=1, color='black')
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(out, dpi=80)


    @staticmethod
    def plotDistributionWithRegressionDinuc(lXs, lYs, rXs, rYs, slope, intercept, R2, pval, meanPhase, stdvPhase,out="", title="", xax="", yax=""):
        """Draw a simple Distribution and add a regression subgraph"""

        fig = plt.Figure(figsize=(20,20))
        fig.suptitle(title, fontsize=32)
        ax = fig.add_subplot(111)
        ax.plot(lXs,lYs)
        annotation_font = {'size':16}
        for i in rYs:
         #   ax.annotate(str(i),xy=(i-10,lYs[i]+10), **annotation_font)
            ax.annotate(str(i+min(lXs)),xy=(i+min(lXs),lYs[i]), **annotation_font)
        ax.scatter([i+min(lXs) for i in rYs],[lYs[i] for i in rYs],color='red',marker='+', s=200)
        axis_font = {'size':'28'}
        ax.set_xlabel(xax, **axis_font)
        ax.set_ylabel(yax, **axis_font)
        ax.tick_params(labelsize=20)
        #ax.set_xlim(0,lXs[-1]+10)
        #ax.set_ylim(min(lYs)-10,max(lYs)+10)
        ax.set_xlim(min(lXs),max(lXs))
        ax.set_ylim(min(lYs),max(lYs))
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



    @staticmethod
    def plotDistributionWithGeneHistogram(lXs, lYs, lZs, lZOthers,out="out.png", title="title", xax="xax", yax="yax", yax2="yax2"):
        """Draw distribution and an histogram for gene coverage"""

        fig = plt.Figure(figsize=(20,20))
        fig.suptitle(title, fontsize=32)
        gs = gridspec.GridSpec(2, 1, width_ratios=[1],height_ratios=[1, 3], hspace=0.1) 
        ax1 = fig.add_subplot(gs[0])
        ax1.plot(lXs,lZOthers)
        lZmin = [0] * len(lZs)
        ax1.vlines(lXs,lZmin,lZOthers, colors='grey', alpha=0.15)
        if max(lZOthers) <= 0:
            ax1.set_ylim(0,1)
        ax1.set_xlim(lXs[0],lXs[-1])
        lZmax = lZs
        lZmin2 = [300] * len(lZs)
        ax2 = fig.add_subplot(gs[1])
        ax2.vlines(lXs,lZmin,lZmax, colors='grey', alpha=0.15)
        ax3 = ax2.twinx()
        ax3.plot(lXs,lYs)
        ax2.set_xlim(lXs[0],lXs[-1])
        ax2.set_ylim(0,max(lZs)+int(max(lZs)*0.05))
        #ax3.set_ylim(min(lYs)-1,max(lYs)+1)
        axis_font = {'size':'28'}
        ax2.set_xlabel(xax, **axis_font)
        ax3.set_ylabel(yax2, **axis_font)
        ax2.set_ylabel(yax, **axis_font)
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(out, dpi=80)


    @staticmethod
    def plotDistributionWithLimits(lXs, llYs, lKClassif,out="out.png", title="title", xax="xax", yax="yax",legend=""):
        """Draw distributions with upper and lower limits"""

        fig = plt.Figure(figsize=(40,20))
        fig.suptitle(title, fontsize=32)
        nbPlots = len(llYs)
        sqrt = int(math.ceil(math.sqrt(nbPlots)))
        ymax = 0.0
        for val in llYs:
            ymax = max(max(val[0]),ymax)
            ymaxCurrent = max(max(val[2]),ymax)
        ymax = ymax*1.05
        xmax = 147
        gs = gridspec.GridSpec(1,2) 
        ax = fig.add_subplot(gs[0])
        gsLimit = gridspec.GridSpecFromSubplotSpec(sqrt,sqrt, subplot_spec=gs[1])
        for i,val in enumerate(llYs):
            ax.plot(lXs,val[0],color=Graphics.lColors[i%25])
            axCurrent = fig.add_subplot(gsLimit[i]) 
            axCurrent.fill_between(lXs, val[1], val[2], alpha=0.35, edgecolor='black', facecolor=Graphics.lColors[i%25])
            axCurrent.set_title("Cluster K{}, (position: {})".format(i,lKClassif[i]))
            axCurrent.fill_between(lXs, val[3], val[4], alpha=0.85, edgecolor='darkgray', facecolor='lightgray')
            axCurrent.plot(lXs,val[0],color=Graphics.lColors[i%25])
            axCurrent.set_ylim(0,ymaxCurrent)
            axCurrent.set_xlim(1,xmax)
            axCurrent.text(10, ymaxCurrent*0.90, "#nucleosomes: {}".format(legend[i]), fontsize=12)
        axis_font = {'size':'28'}
        ax.set_ylim(0,ymax)
        ax.set_xlim(1,xmax)
        ax.legend(["K{}".format(x) for x in range(0,nbPlots)])
        ax.set_title("all nucleosomes", **axis_font)
        ax.set_xlabel(xax, **axis_font)
        ax.set_ylabel(yax, **axis_font)
        ax.tick_params(labelsize=20)
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(out, dpi=80)

    @staticmethod
    def plotDistributionWithLimitsRefine(lXs, llYs, lKClassif,out="out.png", title="title", xax="xax", yax="yax",legend=""):
        """Draw distributions with upper and lower limits"""

        fig = plt.Figure(figsize=(40,20))
        fig.suptitle(title, fontsize=32)
        nbPlots = len(llYs)
        sqrt = int(math.ceil(math.sqrt(nbPlots)))
        ymax = 0.0
        for i,val in enumerate(llYs):
            if lKClassif[i] != "refine":
                ymax = max(max(val[0]),ymax)
                ymaxCurrent = max(max(val[2]),ymax)
        ymax = ymax*1.05
        xmax = 147
        gs = gridspec.GridSpec(1,2) 
        ax = fig.add_subplot(gs[0])
        gsLimit = gridspec.GridSpecFromSubplotSpec(sqrt,sqrt, subplot_spec=gs[1])
        for i,val in enumerate(llYs):
            if lKClassif[i] != "refine":
                ax.plot(lXs,val[0],color=Graphics.lColors[i%25])
                axCurrent = fig.add_subplot(gsLimit[i]) 
                axCurrent.fill_between(lXs, val[1], val[2], alpha=0.35, edgecolor='black', facecolor=Graphics.lColors[i%25])
                axCurrent.set_title("Cluster K{}, (position: {})".format(i,lKClassif[i]))
                axCurrent.fill_between(lXs, val[3], val[4], alpha=0.85, edgecolor='darkgray', facecolor='lightgray')
                axCurrent.plot(lXs,val[0],color=Graphics.lColors[i%25])
                axCurrent.set_ylim(0,ymaxCurrent)
                axCurrent.set_xlim(1,xmax)
                axCurrent.text(10, ymaxCurrent*0.90, "#nucleosomes: {}".format(legend[i]), fontsize=12)
        axis_font = {'size':'28'}
        ax.set_ylim(0,ymax)
        ax.set_xlim(1,xmax)
        ax.legend(["K{}".format(x) for x in range(0,nbPlots)])
        ax.set_title("all nucleosomes", **axis_font)
        ax.set_xlabel(xax, **axis_font)
        ax.set_ylabel(yax, **axis_font)
        ax.tick_params(labelsize=20)
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(out, dpi=80)



    @staticmethod
    def plotDistributionWithLimitsOld(lXs, llYs ,out="out.png", title="title", xax="xax", yax="yax",legend=""):
        """Draw distributions with upper and lower limits"""

        fig = plt.Figure(figsize=(40,20))
        fig.suptitle(title, fontsize=32)
        nbPlots = len(llYs)
        sqrt = int(math.ceil(math.sqrt(nbPlots)))
        ymax = 0.0
        for val in llYs:
            ymax = max(max(val[0]),ymax)
            ymaxCurrent = max(max(val[2]),ymax)
        ymax = ymax*1.05
        xmax = 147
        gs = gridspec.GridSpec(1,2) 
        ax = fig.add_subplot(gs[0])
        gsLimit = gridspec.GridSpecFromSubplotSpec(sqrt,sqrt, subplot_spec=gs[1])
        for i,val in enumerate(llYs):
            ax.plot(lXs,val[0],color=Graphics.lColors[i])
            ax.set_ylim(0,ymax)
            ax.set_xlim(1,xmax)
            axCurrent = fig.add_subplot(gsLimit[i]) 
            axCurrent.fill_between(lXs, val[1], val[2], alpha=0.35, edgecolor='black', facecolor=Graphics.lColors[i])
            for j in val[3:]:
                axCurrent.plot(lXs,j,color='lightgray')
            axCurrent.plot(lXs,val[0],color=Graphics.lColors[i])
            axCurrent.set_ylim(0,ymaxCurrent)
            axCurrent.set_xlim(1,xmax)
#            CV = variation(val[0])
#            axCurrent.text(60, ymaxCurrent*0.7, "CV = {}".format(CV), fontsize=10)
        axis_font = {'size':'28'}
        ax.set_xlabel(xax, **axis_font)
        ax.set_ylabel(yax, **axis_font)
        ax.tick_params(labelsize=20)
        if legend:
            ax.legend(legend)
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(out, dpi=80)


    @staticmethod
    def plotHierarchicalClustering(lXs,distance,labels=[],out="hierarchical_clustering.png", title="", xax="xax", yax="yax"):
        """Draw hierarchical clustering"""

        #fig = plt.Figure(figsize=(20,20))
        fig = plt.Figure()
        #fig.suptitle(title, fontsize=32)
        #ax = fig.add_subplot(111)
        i = dendrogram(
            lXs,
            leaf_rotation=90.,  # rotates the x axis labels
            leaf_font_size=8.,  # font size for the x axis labels
            color_threshold=distance,
            distance_sort='ascending',
            labels=labels
            #distance_sort='descending'
        )
        axis_font = {'size':'28'}
        #ax.set_xlabel(xax, **axis_font)
        #ax.set_ylabel(yax, **axis_font)
        #ax.tick_params(labelsize=20)
        plt.title = title
        plt.xlabel = xax 
        plt.ylabel = yax 
        #plt.show()
        plt.savefig(out,dpi=80, format='png', bbox_inches='tight')
        plt.close()
        #canvas = FigureCanvasAgg(fig)
        #canvas.print_figure(out, dpi=80)


    @staticmethod
    def plotAutocorrelation(lXs, lYs, out="out.png", title="title", xax="xax", yax="yax"):
        """Draw autocorrelation plot, return list of autocorrelation coefficients"""

#        print len(lXs)
#        print len(lYs)
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
        print(len(lRhs))

        fig = plt.Figure(figsize=(20,20))
        fig.suptitle(title, fontsize=32)
        ax = fig.add_subplot(111)
        ax.plot(lXs,lRhs)
        axis_font = {'size':'28'}
        ax.set_xlabel(xax, **axis_font)
        ax.set_ylabel(yax, **axis_font)
        ax.tick_params(labelsize=20)
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(out, dpi=80)

        return lRhs

    @staticmethod
    def plotPValHistogram(lXs, lYs, out="out.png", title="title", xax="xax", yax="yax"):
        """Draw histogram of p-values"""

        fig = plt.Figure(figsize=(20,20))
        fig.suptitle(title, fontsize=32)
        ax = fig.add_subplot(111)
        ax.hist(lXs,lYs)
        axis_font = {'size':'28'}
        ax.set_xlabel(xax, **axis_font)
        ax.set_ylabel(yax, **axis_font)
        ax.tick_params(labelsize=20)
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(out, dpi=80)


