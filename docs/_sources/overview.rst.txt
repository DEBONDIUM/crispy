.. _overview:

Overview
===============

About
-----

**CRISPY** is a Python package for detecting, visualizing and analyzing fragmentation in bonded particle based simulations (SPH, DEM, Peridynamic).

It provides tools to:

- Read simulation files.
- Detect connected fragments (2D and 3D).
- Visualize the evolution of fragments over time.
- Filter out small trash fragments and manage fragmentation hierarchy.
- Analyze fragment properties (sphericity, volume, area, ...)

Main Features
-------------

- 📁 Load simulation data from bonded particle based simulation files.
- 🧩 Detect and extract fragments based on node connectivity.
- 🕸️ Build and analyze the hierarchy of parent-child fragments.
- 🧼 Automatically discard small trash fragments using a threshold.
- 📊 Interactive 2D/3D visualization using PyVista.
- 🧪 Easily integrated into analysis pipelines or used as a standalone tool.

Objectives
----------

- Provide a robust framework to extract and analyze fragmentation patterns (2D and 3D).
- Facilitate the visualization of fragment evolution across iterations.
- Offer an intuitive Python interface for post-simulation analysis.
- Allow researchers to explore fragment distribution over severe loading condition simulation.

