
ZEBRAFISH
===============================================================================

Please see ../eltopo3d/readme.txt for information on building the El Topo 
library.

This is an application based on El Topo library to simulate invagination.  This readme describes:

- how to build the executable
- how to run the executable
- an overview of the executable source code

Note: This is research-grade code -- full of hacks, bugs, dead code, etc.

Building Zebrafish:
=====================

1) Create Makefile.local_defs

The included Makefile reads a file called Makefile.local_defs, which contains 
platform-specific definitions.  You must create this file!  
Makefile.example_defs includes suggested settings for Linux and OS/X.

If you want to use the command-line interface, define NO_GUI in your build 
(e.g. in Makefile.local_defs).  If you want to use the GUI, be sure to link 
against OpenGL and GLUT.  

2) Generate dependencies and compile

Building the executable is done by running "make depend" followed by 
"make release" or "make debug".  This should also automatically build the
El Topo library.

Example:
$> make depend
$> make release

Using Zebrafish:
=====================

Launching the executable requires two command line parameters:
1) A path to a text "script" file specifying the simulation type and initial 
geometry
2) A path specifying where output files will be written

Example:
$> ./zebrafish_release scripts/simulation-parameters.txt /var/tmp/

If you are running the GUI, this should pop up a GLUT window with a view of the 
triangle mesh surface (see main.cpp to learn what the keyboard does).  If you 
are running with no GUI, it will immediately start running the simulation 
defined by the script, outputting one binary mesh file per frame in /var/tmp/.


Code base:
=====================

Summaries of the important source files used by El Topo are found in the readme
in that directory.

Drivers:
---------------------

In zebrafish, the dynamic surface is assigned a linear velocity per vertex by a 
"mesh driver".

References
=====================

[Bridson et al. 2007]: R. Bridson, J. Hourihan, and M. Nordenstam,  Curl noise 
for procedural fluid flow, Proc. ACM SIGGRAPH 2007

[Brochu and Bridson 2009]: Tyson Brochu and Robert Bridson, Robust Topological 
Operations for Dynamic Explicit Surfaces, SIAM J. Sci. Comput., vol. 31, no. 4 
(2009), pp. 2472-2493 

[Enright et al. 2002]: D. Enright, R. Fedkiw, J. Ferziger, and I. Mitchell, A 
hybrid particle level set method for improved interface capturing, J. Comput. 
Phys., 183 (2002), pp. 83???116.

[Jiao 2007]: X. Jiao, Face offsetting: A unified framework for explicit moving 
interfaces, J. Comput. Phys., 220 (2007), pp. 612???625.


