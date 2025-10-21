#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 14:26:21 2025

@author: lbremaud
"""

# import crispy
import crispy as cp

# load files
files = cp.load_files("files/", extension="agdd")

# detect fragments from stacked files
frag_dict = cp.detect_fragment(files)

# set camera position
camera_position=[
    (-3.02e-05, -9.30e-05, 0.14),
    (-3.02e-05, -9.30e-05, 0.0),
    (0.0, 1.0, 0.0)]

# plot fragments for all iterations
for it in frag_dict.keys():
    cp.pvplot(frag_dict[it], filename=f"plot2D_it{it}", auto_close=True, camera_position=camera_position)

# stackplot of the fragments over iterations
cp.stackplot([f for fgroup in frag_dict.values() for f in fgroup])

# statistic histogram of fragments at iteration n°50
cp.histplot(frag_dict[49])

# statistics of fragments at iteration n°50
cp.stats(frag_dict[49])



