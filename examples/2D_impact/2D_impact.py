#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 14:26:21 2025

@author: lbremaud
"""

# =============================================================================
# LIB
# =============================================================================
import crispy as cp

# =============================================================================
# MAIN
# =============================================================================
# Initialize FragmentDetector and specified .aggd files path
detector = cp.FragmentDetector("agdd/")

# Build fragments
detector.build_fragments()

# 2D plot at each iteration
for it in range(detector.iterations_nb):
    detector.plot2D(iteration=it, save=True, save_format="eps")

# Stack plot: repartition of area
detector.stackplot()

# Graph: heredity
detector.graphplot()