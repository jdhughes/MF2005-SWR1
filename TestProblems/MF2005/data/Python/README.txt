09/12/2012

The python scripts in this subdirectory can be used to extract data
from the SWR1 output files for Test Simulations 1-5 documented in the
SWR1 Process report (http://pubs.usgs.gov/tm/6a40/).

Included scripts:
    SWRSample01.py   
    SWRSample02.02.py
    SWRSample02.py   
    SWRSample03.py   
    SWRSample04.py   
    SWRSample05.py   

Scripts can be run from the command line by typing:
    python SWRSample01.py 

Alternatively, the provided batch file (PlotAll.bat) can be executed.

Output from the scripts is written to the ..\Figures\ subdirectory.

The scripts have been tested using python2.6 but will also work with
python2.7. Unknown if scripts will work with python3.0.

Python dependencies include:
    NUMPY
    MATPLOTLIB
    
Scripts also depend on included MFBinaryClass.py. This script includes
functionality for reading all binary SWR1 output for version 1.0 of
the SWR1 Process.

Only limited assistance can be provided for running the provided scripts
or applying MFBinaryClass.py to extract binary SWR1 Process output for
specific problems.

Point of contact:
    Joseph D. Hughes
    Florida Water Science Center
    4446 Pet Lane
    Suite 108
    Lutz, FL 33559
    jdhughes@usgs.gov 