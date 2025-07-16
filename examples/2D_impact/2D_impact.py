#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 14:26:21 2025

@author: lbremaud
"""

import crispy as cp

sim = cp.FragmentDetector("agdd/")
sim.build_fragments()

# for i in range(sim.iterations_nb):
    # sim.plot2D(iteration=i, save = True)
sim.viewer3D()