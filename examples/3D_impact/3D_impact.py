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
    (-0.25, 0.05, 0.21),
    (0.0, 0.0, 0.0),
    (0.19, 0.98, 0.0)]

# plot fragments for all iterations
for it in frag_dict.keys():
    cp.pvplot(frag_dict[it], filename=f"plot3D_it{it}", camera_position=camera_position)





