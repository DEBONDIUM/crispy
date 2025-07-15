#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 14:26:21 2025

@author: lbremaud
"""
# import sys
# from pathlib import Path
# sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import crispy as cp

sim = cp.Crispy("agdd/")
sim.build_fragments()

for i in range(sim.iterations_nb):
    sim.plot2D(iteration=i, save = True)
sim.viewer3D()