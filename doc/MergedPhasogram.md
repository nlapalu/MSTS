Here a simple script to generate a merged phasogram from several runs of MSTS_feature_phasogram.py 

write in list.txt all the *.phaso file to merge with the associated legend (tab delimited).

__*list.txt*__:

```
1.phaso	1TPM>X
5-1.phaso	5TPM>X>1TPM
50-5.phaso	50TPM>X>5TPM
50.phaso	X>50TPM
```

write the code below in mergePhasograms.py (up to 6 conditions) and run the command:

`python mergePhasograms.py list.txt`

Your phasogram will be exported in phaso.merge.png

__*script mergePhasograms.py:*__

```python
#!/usr/bin/env python

import sys

import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg

color=["red","blue","green","orange","purple","maroon"]
llXs = []
llYs = []
legend = []
with open(sys.argv[1]) as f:
    for line in f:
        values = line.rstrip().split('\t')
        legend.append(values[1])
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
                print "Error indice list"
                sys.exit(1) 
        else:
            llXs = lXs


fig = plt.Figure(figsize=(20,20))
fig.suptitle("merged phasogram", fontsize=32)
ax = fig.add_subplot(111)
for i,val in enumerate(llYs):
    ax.plot(lXs,val,color=color[i])
axis_font = {'size':'28'}
ax.set_xlabel("window, bp", **axis_font)
ax.set_ylabel("signal coverage", **axis_font)
ax.tick_params(labelsize=20)
ax.legend(legend, fontsize=20)
canvas = FigureCanvasAgg(fig)
canvas.print_figure("phaso.merge.png", dpi=80)
```
