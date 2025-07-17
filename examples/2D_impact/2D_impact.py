#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 14:26:21 2025

@author: lbremaud
"""

import crispy as cp

detector = cp.FragmentDetector("agdd/")
detector.build_fragments()

for it in range(detector.iterations_nb):
    detector.plot2D(iteration=it)

detector.stackplot()
detector.graphplot()