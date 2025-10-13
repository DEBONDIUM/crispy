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
    def centroid(self) -> Node:
        xyz = np.array([(n.x, n.y, n.z) for n in self._nodes])
        cx, cy, cz = xyz.mean(axis = 0)
        return Node(-1, cx, cy, cz)



















