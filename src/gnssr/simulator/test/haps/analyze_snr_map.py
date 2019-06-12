#!/usr/bin/env python

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import gnssr.simulator.rcs.target_rcs as target_rcs
import gnssr.simulator.rcs.sea_rcs as sea_rcs
from gnssr.simulator.isolines import *
from gnssr.simulator.simulation_configuration import *
from gnssr.simulator.ddm import *

from gnssr.utils import *
from gnssr.simulator.jacobian.planar import *
#!/usr/bin/env python

import pickle

# Getting back the objects:
with open('results_snr_map.pkl', 'rb') as f:
    x_grid, y_grid, z_grid = pickle.load(f)

    np.place(z_grid, z_grid < -30, -30haps_20ms_snr_mapo)

    fig_snr_grid, ax_snr_grid = plt.subplots(1,figsize=(10, 4))
    #contour = plt.contourf(x_grid, y_grid, z_grid, cmap="jet")
    contour = ax_snr_grid.imshow(z_grid, cmap="jet",
            origin="lower",
            extent=(
                x_grid[0][0], x_grid[0][-1],
                y_grid[0][0], y_grid[-1][0],
                ),
            aspect = "auto"
            )
    plt.xlabel("[km]")
    plt.ylabel("[km]")
    fig_snr_grid.colorbar(contour, label="SNR [dB]")

    ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/1000))
    ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000))
    ax_snr_grid.xaxis.set_major_formatter(ticks_x)
    ax_snr_grid.yaxis.set_major_formatter(ticks_y)

    plt.show()
