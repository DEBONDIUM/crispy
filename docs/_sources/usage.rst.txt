.. _usage:

Quick Start
===========

Getting Started
---------------

To begin using Crispy, import the library:

.. code-block:: python

   import crispy as cp

Then load input files from a specified directory:

.. code-block:: python

   files = cp.load_files("files/")

Then launch the fragment detection:

.. code-block:: python

   frag_dict = cp.detect_fragment(files)

Input File Format
-----------------

Crispy operates on all type of files, as long as they contain structured data describing bonded particle based domain properties such as nodes and bonds.

A typical input file should be structured like this:

.. code-block:: text

   # Example input file
   1000 # Number of nodes (Pos X / Pos Y / Pos Z)
   0.	0.	0.
   1.	1.	1.
   2.	2.	2.
   ...
   5000 # Number of bonds (Node 1 / Node 2)
   0	1
   1	2
   3	2
   ...

Make sure all input files in your directory follow a consistent structure.


