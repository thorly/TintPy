#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import colorcet as cc
import re

#from colorcet.plotting import swatch, swatches
#import holoviews as hv
#hv.extension('matplotlib')
#sws = swatches()
#hv.save(sws, 'swatches.png', fmt='png', size=1000)

cgrps=['glasbey','diverging','linear', 'cyclic', 'isoluminant','rainbow']

for cgrp in cgrps:
  cmap_list=cc.all_original_names(group=cgrp)
  print('\n%s:'%cgrp)
  for cmn in cmap_list:
    cmlist = re.split(', ', cc.get_aliases(cmn))
    print(cmn, ':', cmlist)
    for c0 in cmlist:
      c1=c0.strip()
      print('  colormap: ',c1)
      cm = cc.cm[c1]
      cmf = open(c1 + '.cm', 'w')
      for j in range(256):
        cl=cm(j)
        r=int(round(cl[0]*255.0))
        g=int(round(cl[1]*255.0))
        b=int(round(cl[2]*255.0))
        cmf.write("%3d %3d %3d\n"%(r,g,b))
      cmf.close()
