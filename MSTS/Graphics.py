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


