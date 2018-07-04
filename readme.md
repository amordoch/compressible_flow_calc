Compressible Flow Calculator (CFC)

Copyright (C) 2018 Ariel Mordoch

See LICENSE for license details
## What is this?
This Python program was built to solve problems related to compressible flow aerodynamics, specifically: isentropic flow
with area change, normal and oblique shock analysis, and Prandtl-Meyer expansion fans. 

Several other libraries do already exist which do the same thing as
CFC, but better [(caeroc)](https://github.com/ashwinvis/caeroc). 
CFC is mostly intended to be a personal project and 
learning experience.
## Installation
Installation is pretty simple. Download either the source code from the releases section
or the wheel/tarball directly out of `dist/` and run the following
in your terminal:
- If you downloaded the wheel:

  `pip install compressible_flow_calc-0.1-py3-non-any.whl` 
 - If you downloaded the tarball:
 
   `pip install compressible_flow_calc-0.1.tar.gz`
 
 A release to PyPi is planned.
 ## Requirements
 CFC requires numpy. If you install with pip it will be installed
 automatically if you don't already have it. 
## Usage
You can use compressible_flow_calc in several ways:
1. Directly from the command line: 
   
   `python -m compressible_flow_calc`
2. As a package:

   >yourmodule.py
   ```python
   >>> import compressible_flow_calc as cfc
   >>> cfc.calc.A_over_Astar(M=2, gamma=1.4)
   1.6875
   ```
3. Directly executing cli.py:

   `you@yourcomputer:/path/to/lib$ python cli.py`
 
 For options 1 and 3, CFC has its own "CLI" built in which will
 allow you to do simple mathematical analyses. For more complex programs,
 CFC exposes its mathematical functions as an importable Python package. 
 ## Documentation
 ...Currently doesn't really exist. But so far, there is not much functionality
 anyway so feel free to take a look at the comments.
 ## A Note on Units
 So far, all functions implemented are ratios -- which are unitless
 by definition. Use whatever unit system you prefer. 
 
 Theorectically there should be no problems with using values with
 units defines by pint, but this is not tested.