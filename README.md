#### alps_cthyb  
`alps_cthyb` is the hybridization expansion continuous-time quantum Monte Carlo code, implemented and described in: 
H. Hafermann, P. Werner, E. Gull, CPC 184, 1280 (2013).
This version of the code is its adaptation to *ALPSCore* library. 

#### Dependencies
1. ALPSCore (http://alpscore.org)
2. LAPACK

#### Usage 
``
alps_cthyb parameters
``
See `alps_cthyb --help` for the list of parameters

#### Interface changes from the ALPS version 
Parameters:
- MAX_TIME -> timelimit

Python bindings are temporarily disabled

