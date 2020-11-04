# minkmaps
**C++ code for a Morphometric Statistical Analysis via Minkowski Maps**  
Author: Michael Andreas Klatt (<software@mklatt.org>)  
License: GNU GPLv3  
Version: 19.06  

The morphology of a binned map of events is quantified by Minkowski
functionals and a null hypothesis test is performed whether the events
are completely random (independent and uniformly distributed).

Due to its first application, a simulated or measured count map is also
called sky map.

**Based on the paper**  
Michael Andreas Klatt, Klaus Mecke.
[Detecting structured sources in noisy images via Minkowski maps](https://dx.doi.org/10.1209/0295-5075/128/60001), EPL (Europhysics Letters) 128:60001 (2019)

Please cite the paper if you use the code.

Compilation
-----------
Via makefile, tested on linux, requires Boost and GSL.

Usage
-----
**minkmap.cpp**
contains main function and is used to set up the simulation of a sky map
or the reading or input data, as well as the statistical analysis.

**minkmap.conf**
is a configuration file (can be renamed) to store basic parameters.
Alternatively command line options can be used, for example:

                       ./minkmap -l 5 -U 75 -V 33 -w 15 -k 7 -s 17

Known issues
------------

* Including header-files etc. contains "historic constructions" and
  should be improved.
* A proper documentation is missing. If required (software@mklatt.org),
  this can easily be provided by turning the comments into doxygen
  commands.
* Option for minus sampling boundary conditions (mbc) in binnedskymap is deprecated

