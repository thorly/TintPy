#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmocean

cmaps = cmocean.cm.cmap_d
for cmn in cmaps.keys():
  cm=cmaps[cmn]
  print(cmn+'.cm')
  cmf = open(cmn + '.cm', 'w')
  for j in range(256):
    cl=cm(j)
    r=int(round(cl[0]*255.0))
    g=int(round(cl[1]*255.0))
    b=int(round(cl[2]*255.0))
    cmf.write("%3d %3d %3d\n"%(r,g,b))
  cmf.close()
