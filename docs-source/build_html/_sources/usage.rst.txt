.. _usage:

Quick Start
===========

Getting Started
---------------

To begin using Crispy, import the library and initialize a fragment detector:

.. code-block:: python

   import crispy as cp

   detector = cp.FragmentDetector("path/to/files")
   detector.build_fragments()


AGDD File Format
----------------

Crispy operates on `.agdd` files, which contain structured data describing discrete domain properties such as nodes and bonds.

A typical `.agdd` file looks like this:

.. code-block:: text

   # Example AGDD file
   1000 # Number of nodes (Pos X / Pos Y / Pos Z / Radius)
   0.	0.	0.	0.5
   1.	1.	1.	0.6
   2.	2.	2.	0.4
   ...
   5000 # Number of bonds (Node 1 / Node 2)
   0	1
   1	2
   3	2
   ...

Make sure all `.agdd` files in your directory follow a consistent structure.


