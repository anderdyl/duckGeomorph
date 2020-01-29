#!/usr/bin/env python

# Copyright 2013 Perttu Luukko

# This file is part of libeemd.
# 
# libeemd is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# libeemd is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with libeemd.  If not, see <http://www.gnu.org/licenses/>.

from pyeemd import ceemdan
from pyeemd.utils import plot_imfs
from matplotlib.pyplot import plot, show, title
import numpy as np

# Decompose a delta-function signal
N = 512
signal = np.zeros(N)
signal[N//2] = 1

# Plot the original data using Matplotlib
title("Original signal")
plot(signal)

# Calculate IMFs and the residual by CEEMDAN
imfs = ceemdan(signal, noise_strength=0.2, ensemble_size=500)

# Plot the results using the plot_imfs helper function from pyeemd.utils
plot_imfs(imfs, plot_splines=False)
show()

# You can compare the results with Fig. 1 in the original CEEMDAN paper at
# http://dx.doi.org/10.1109/ICASSP.2011.5947265
