#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

xx=np.arange(0,10,0.01)
yy=xx*np.exp(-xx)

path = Path(np.array([xx,yy]).transpose())
patch = PathPatch(path, facecolor='none')
plt.gca().add_patch(patch)

im = plt.imshow(xx.reshape(yy.size,1),  cmap='jet',
                origin='lower',extent=[0,10,-0.0,0.40],aspect="auto", clip_path=patch, clip_on=True)

plt.show()
