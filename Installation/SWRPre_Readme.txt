

                       SWRPre - Version: 1.0.0
                  A Preprocessor for the SWR1 Process


NOTE: Any use of trade, product or firm names is for descriptive purposes 
      only and does not imply endorsement by the U.S. Government.

SWRPre version 1.0 is packaged for personal computers using one of the 
Microsoft Windows operating systems. An executable file for personal 
computers are provided. The executable files were compiled on a personal
computer with the Intel(R) Core(TM) I7-2820QM CPU chipset, running the 
Microsoft Windows 7 Enterprise operating system, using the Microsoft 
Visual Studio 2008 Version 9.0.21022.8 development environment and the 
Microsoft Visual C# 2008 compiler. 

Instructions for installation, execution, and testing of this version of
SWRPre are provided below.



                            TABLE OF CONTENTS

                         A. DISTRIBUTION FILE
                         B. INSTALLING
                         C. EXECUTING THE SOFTWARE
                         D. TESTING
                         E. COMPILING


A. DISTRIBUTION FILE

The following self-extracting distribution file is for use on personal
computers:

         SWRPre_Installation.exe

The distribution file contains:

          Executable for SWRPre.
          MF2005-SWR documentation.

The distribution file is a self-extracting program.  Execution of the
distribution file creates numerous individual files.  The extraction
program allows you to specify the directory in which the files should
be restored.  The default installation directory is C:\WRDAPP.  You
have the opportunity to specify an alternate installation directory
during extraction of the software. The following directory structure
will be created in the installation directory:


   |
   |--SWRPre
   |    |--bin           ; Compiled SWRPre executable for personal computers
   |    |--doc           ; Documentation report for MF2005-SWR

Included in directory SWRPre\doc is the SWR1 Process documentation 
report, which is a Portable Document Format (PDF) file. The PDF file is 
readable and printable on various computer platforms using Acrobat Reader 
from Adobe. The Acrobat Reader is freely available from the following World
Wide Web site:

      http://www.adobe.com/


B. INSTALLING

To make the executable version of SWRPre accessible from any
directory, the directory containing the executables (SWRPre\bin)
should be included in the PATH environment variable (see explanation below).  

As an alternative, the executable file, SWRPre.exe and MF2005-SWR_x64.exe, 
in the SWRPre\bin directory can be copied into a directory already
included in the PATH environment variable.

       HOW TO ADD TO THE PATH ENVIRONMENT VARIABLE
          WINDOWS 9 X AND WINDOWS ME SYSTEMS
          
Add the following line to the AUTOEXEC.BAT file:

  PATH=%PATH%;C:\WRDAPP\SWRPre\bin

Note, reboot your system after modifying AUTOEXEC.BAT.


       HOW TO ADD TO THE PATH ENVIRONMENT VARIABLE
               WINDOWS NT SYSTEMS

From the Start menu, select Settings and then Control Panel.  Double click
System and select the Environment tab. To add a new user variable, enter
"PATH" in the Variable field and enter

   %PATH%;C:\WRDAPP\SWRPre\bin

in the Value field.  Click Set and then click OK.  If a PATH user variable
already is defined, click on it in the User Variables pane, add
";C:\WRDAPP\SWRPre\bin" to its definition in the Value field, and click
OK.  Initiate and use a new Windows Command Prompt window after making this
change.


       HOW TO ADD TO THE PATH ENVIRONMENT VARIABLE
             WINDOWS 2000 OR XP SYSTEMS
             
From the Start menu, select Settings and then Control Panel.  Double click
System and select the Advanced tab.  Click on Environment Variables.  If
a PATH user variable already is defined, click on it in the User Variables
pane, then click Edit.  In the Edit User Variable window, add
";C:\WRDAPP\SWRPre\bin" to the end of the Variable Value (ensure that
the current contents of the User Value are not deleted) and click OK.  If
a PATH user variable is not already defined, in the User variables pane of
the Environment Variables window, click New.  In the New User Variable
window, define a new variable PATH as shown above.  Click OK.  Click OK
in the Environment Variables window and again in the System Properties
window.  Initiate and use a new Windows Command Prompt window.


C. EXECUTING THE SOFTWARE

After the executable files in the SWRPre\bin directory are installed in
a directory that is included in your PATH, can be executed from the 
Windows Command-Prompt window using the commands:

      SWRPre.exe

SWRPre can also be excuted by double-clicking on the executable. See the 
documenation (TM6-A40) for more information on use of SWRPre. 


D. TESTING

No sample problems are available.


E. COMPILING

The executable files provided in SWRPre\bin were created using the Microsoft 
Visual C# 2008 compiler. Although an executable version of the program is  
provided, the source code can be provided so that SWRPre can be recompiled 
if necessary. However, the USGS cannot provide assistance to those compiling 
SWRPre.
