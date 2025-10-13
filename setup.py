#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 14:47:08 2025

@author: lbremaud
"""

from setuptools import setup, find_packages

setup(
    name="crispy",
    version="0.2.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "matplotlib",
        "scipy",
        "open3d",
        "networkx",
    ],
    author="DEBONDIUM",
    description="Fragment detection and visualization from AGDD files",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/DEBONDIUM/crispy",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)
