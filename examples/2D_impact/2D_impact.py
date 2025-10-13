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
fragments_history = cp.detect_fragment(files)

# stackplot
cp.stackplot(fragments_history, save=False)

# plot fragments for all iterations
for it in fragments_history.keys():
    cp.plot(fragments_history[it], filename=f"plot2D_it{it}")
