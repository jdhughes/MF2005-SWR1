

                      MF2005-SWR - Version: 1.01.0
                  MODFLOW-2005 with the SWR1 Process


NOTE: Any use of trade, product or firm names is for descriptive purposes 
      only and does not imply endorsement by the U.S. Government.

MF2005-SWR version 1.0 is packaged for personal computers using one of the 
Microsoft Windows operating systems. Executable files for personal 
computers are provided as well as the source code. The executable files were 
compiled on a personal computer with the Intel(R) Core(TM) I7-2820QM CPU 
chipset, running the Microsoft Windows 7 Enterprise operating system, 
using the Microsoft Visual Studio 2008 Version 9.0.21022.8 development
environment and the Intel(R) Visual Fortran Version 12.1.2.278 compiler. 
The source code is provided to aid users in compilation on 
other computers. However, no support is provided for compilation.

IMPORTANT: Users should review the file Summary_MF2005-SWR.txt for a description
of, and references for, this software. Users should also review the file 
release.txt, which describes changes that have been introduced into MF2005-SWR
with each official release; these changes may substantially affect users.

Instructions for installation, execution, and testing of this version of
MF2005-SWR are provided below.



                            TABLE OF CONTENTS

                         A. DISTRIBUTION FILE
                         B. INSTALLING
                         C. EXECUTING THE SOFTWARE
                         D. TESTING
                         E. COMPILING


A. DISTRIBUTION FILE

The following self-extracting distribution file is for use on personal
computers:

         mf2005-SWR_Installation.exe

The distribution file contains:

          Executables and source code for MF2005-SWR.
          Executable for SWRPre.
          MF2005-SWR documentation.
          Five MF2005-SWR sample problems.

The distribution file is a self-extracting program.  Execution of the
distribution file creates numerous individual files.  The extraction
program allows you to specify the directory in which the files should
be restored.  The default installation directory is C:\WRDAPP.  You
have the opportunity to specify an alternate installation directory
during extraction of the software. The following directory structure
will be created in the installation directory:


   |
   |--MF2005-SWR_1.9
   |    |--bin              ; Compiled MF2005-SWR executables for personal computers
   |    |--data             ; Sample problems
   |    |--doc              ; Documentation report for MF2005-SWR
   |    |--output_test      ; Output files from running the sample problems
   |    |--output_test_x64  ; Output files from running the sample problems with
   |                          64-bit executable
   |    |--src              ; Source code for MF2005-SWR

Included in directory MF2005-SWR_1.9\doc is the SWR1 Process documentation 
report, which is a Portable Document Format (PDF) file. The PDF file is 
readable and printable on various computer platforms using Acrobat Reader 
from Adobe. The Acrobat Reader is freely available from the following World
Wide Web site:

      http://www.adobe.com/


B. INSTALLING

To make the executable versions of MF2005-SWR accessible from any
directory, the directory containing the executables (MF2005-SWR_1.9\bin)
should be included in the PATH environment variable (see explanation below).  

As an alternative, the executable files, MF2005-SWR.exe and MF2005-SWR_x64.exe, 
in the MF2005-SWR_1.9\bin directory can be copied into a directory already
included in the PATH environment variable.

       HOW TO ADD TO THE PATH ENVIRONMENT VARIABLE
          WINDOWS 9 X AND WINDOWS ME SYSTEMS
          
Add the following line to the AUTOEXEC.BAT file:

  PATH=%PATH%;C:\WRDAPP\MF2005-SWR_1.9\bin

Note, reboot your system after modifying AUTOEXEC.BAT.


       HOW TO ADD TO THE PATH ENVIRONMENT VARIABLE
               WINDOWS NT SYSTEMS

From the Start menu, select Settings and then Control Panel.  Double click
System and select the Environment tab. To add a new user variable, enter
"PATH" in the Variable field and enter

   %PATH%;C:\WRDAPP\MF2005-SWR_1.9\bin

in the Value field.  Click Set and then click OK.  If a PATH user variable
already is defined, click on it in the User Variables pane, add
";C:\WRDAPP\MF2005-SWR_1.9\bin" to its definition in the Value field, and click
OK.  Initiate and use a new Windows Command Prompt window after making this
change.


       HOW TO ADD TO THE PATH ENVIRONMENT VARIABLE
             WINDOWS 2000 OR XP SYSTEMS
             
From the Start menu, select Settings and then Control Panel.  Double click
System and select the Advanced tab.  Click on Environment Variables.  If
a PATH user variable already is defined, click on it in the User Variables
pane, then click Edit.  In the Edit User Variable window, add
";C:\WRDAPP\MF2005-SWR_1.9\bin" to the end of the Variable Value (ensure that
the current contents of the User Value are not deleted) and click OK.  If
a PATH user variable is not already defined, in the User variables pane of
the Environment Variables window, click New.  In the New User Variable
window, define a new variable PATH as shown above.  Click OK.  Click OK
in the Environment Variables window and again in the System Properties
window.  Initiate and use a new Windows Command Prompt window.


C. EXECUTING THE SOFTWARE

A 32 bit and a 64 bit executable are provided in the MF2005-SWR_1.9\bin  
directory. Two executables are provided because computers often use either  
the 32 bit Windows XP or the 64 bit Windows 7 operating systems. Large 
simulations may not run on a 32 bit operating system due to limitations 
in the amount of available random access memory (RAM). A 64 bit operating
system provides much more available RAM than a 32 bit operating system. 
Thus, it is recommended that a 64 bit executable be used on a 64 bit operating  
system for large simulations.   

After the executable files in the MF2005-SWR_1.9\bin directory are installed in
a directory that is included in your PATH, MF2005-SWR is initiated in
a Windows Command-Prompt window using the commands:

      mf2005-SWR.exe [Fname]

or
      mf2005-SWR_x64.exe [Fname]

The optional Fname argument is the name of the MF2005-SWR Name File.

The data arrays in MF2005-SWR are dynamically allocated, so models are not
limited by hard-coded array limits. However, it is best to have at least 
2 MB of RAM available to hold all of the required data. If there is 
less available RAM than the model requires, which depends on the size of the 
application, the program will use virtual memory; however, this can
slow execution significantly. If there is insufficient RAM to run
the model, then MF2005-SWR will not initiate the beginning of the 
simulation; however, the Windows Command-Prompt window may continue to 
indicate that MF2005-SWR is executing. For this circumstance, the program 
must be terminated manually using the Windows Task Manager application.

Some of the files written by MF2005-SWR are unformatted files. The structure
of these files depends on the compiler and options in the Fortran write
statement.  MF2005-SWR is compiled with the unformatted file type specified
as "BINARY". Any program that reads the unformatted files produced by
MF2005-SWR must be compiled with a compiler that produces programs that
use the same structure for unformatted files.  For example, Zonebudget
and Modpath use unformatted budget files produced by MF2005-SWR.  
Another example is head files that are generated by one MF2005-SWR
simulation and used in a following simulation as initial heads.  Both 
simulations must be run using an executable version of MF2005-SWR that uses 
the same unformatted file structure.

A simple preprocessor (SWRPre) for creating SWR1 Process connectivity 
from an Esri polyline shapefile and a MODFLOW discretization (DIS) file
is also included. SWRPre creates a data file that defines SWR1 reaches, 
reach groups, and reach connectivity. SWRPre also creates (1) an Esri 
polyline shapefile that contains the data contained in the SWR1 data file 
to validate the SWR1 network, and (2) an Esri polygon shapefile of the 
DIS file to validate the coordinate offset and rotation angle provided 
to the pre-processor. SWRPre is installed in the MF2005-SWR_1.9\bin 
directory and can be executed from the Windows Command-Prompt window 
using the commands:

      SWRPre.exe

SWRPre can also be excuted by double-clicking on the executable. See the 
documenation (TM6-A40) for more information on use of SWRPre. 


D. TESTING

Five sample problems with MF2005-SWR_1.9 data sets are provided to verify that  
MF2005-SWR_1.9 is correctly installed and running on the system.  The sample 
problems also may be looked at as examples of how to use the program.


E. COMPILING

The executable files provided in MF2005-SWR_1.9\bin were created using the Intel  
Visual Fortran 12.1.2.278 compiler. Although executable versions of the program are  
provided, the source code also is provided in the MF2005-SWR_1.9\src directory so 
that MF2005-SWR_1.9 can be recompiled if necessary. However, the USGS cannot provide  
assistance to those compiling MF2005-SWR_1.9.
