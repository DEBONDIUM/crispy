#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  7 14:57:21 2025

@author: lbremaud
"""

# =============================================================================
# LIB
# =============================================================================
import numpy as np
from scipy.spatial import ConvexHull

# =============================================================================
# CLASS NODE
# =============================================================================
class Node:
    def __init__(self, i: int, x: float, y: float, z: float):
        self._i = int(i)
        self._x = float(x)
        self._y = float(y)
        self._z = float(z)
    
    # ---- Properties ----
    @property
    def i(self) -> int:
        return self._i

    @property
    def x(self) -> float:
        return self._x
    
    @property
    def y(self) -> float:
        return self._y
    
    @property
    def z(self) -> float:
        return self._z

# =============================================================================
# CLASS BOND
# =============================================================================
class Bond:
    def __init__(self, i: int, n1: Node, n2: Node):
        self._i = int(i)
        self._n1 = n1
        self._n2 = n2
    
    # ---- Properties ----
    @property
    def i(self) -> int:
        return self._i

    @property
    def n1(self) -> Node:
        return self._n1
    
    @property
    def n2(self) -> Node:
        return self._n2

# =============================================================================
# CLASS FRAGMENT
# =============================================================================
class Fragment:
    def __init__(self, i: int, nodes: list[Node], bonds: list[Bond], it: int,
                 ratio: float, ancestor: list[int] = []):
        self._i = int(i)
        self._nodes = nodes
        self._bonds = bonds
        self._it = int(it)
        self._ratio = float(ratio)
        self._ancestor = ancestor
        
        tol = 1e-8
        xyz = np.array([(n.x, n.y, n.z) for n in self._nodes])
        dim_mask = np.ptp(xyz, axis=0) > tol
        pts = xyz[:, dim_mask]

        cx, cy, cz = xyz.mean(axis = 0)
        self._centroid = Node(-1, cx, cy, cz)
        
        self._dim = np.sum(dim_mask)
        if self._dim == 3:
            hull = ConvexHull(pts)
            self._volume = hull.volume
            self._area = hull.area
            self._perimeter = np.nan
            self._sphericity = (np.pi ** (1/3) * (6 * self._volume) ** (2/3)) / self._area
        elif self._dim == 2:
            hull = ConvexHull(pts)
            self._volume = np.nan
            self._area = hull.volume
            self._perimeter = hull.area
            self._sphericity = 4 * np.pi * self._area / (self._perimeter ** 2)
        else:
            self._volume = np.nan
            self._area = np.nan
            self._perimeter = np.nan
            self._sphericity = np.nan
    
    # ---- Properties ----
    @property
    def i(self) -> int:
        return self._i
    
    @property
    def nodes(self) -> list[Node]:
        return self._nodes
    
    @property
    def bonds(self) -> list[Bond]:
        return self._bonds
    
    @property
    def it(self) -> int:
        return self._it
    
    @property
    def ratio(self) -> float:
        return self._ratio
    
    @property
    def ancestor(self) -> list[int]:
        return self._ancestor
    
    @property
    def dim(self) -> int:
        return self._dim
    
    @property
    def volume(self) -> float:
        return self._volume
    
    @property
    def area(self) -> float:
        return self._area
    
    @property
    def perimeter(self) -> float:
        return self._perimeter
    
    @property
    def sphericity(self) -> float:
        return self._sphericity
    
    @property
    def centroid(self) -> Node:
        return self._centroid



















