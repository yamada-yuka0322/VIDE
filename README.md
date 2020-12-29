```
\        /   /   |-\    -----
 \      /    |   |  \   |
  \    /    /    |   |  |--
   \  /     |    |  /   |
    \/      /    |-/    -----

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
```

This is VIDE, the Void IDentification and Examination toolkit.

For more information, see http://www.cosmicvoids.net

Please cite arXiv:1406.1191 and arXiv:0712.0349 if you use this software, 
using the following suggested sentence:

"This work uses voids identified with VIDE\footnote{\url{http:www.cosmicvoids.net}} (Sutter et al. 2014), 
which implements an enhanced version of ZOBOV (Neyrinck 2008) to construct 
voids with a watershed algorithm."


License/Copyright information
-----------------------------

Copyright (C) 2010-2020 Guilhem Lavaux, 2011-2014 P.M. Sutter.
This software is put under the GNU Public License. 
Please see LICENSE for further information.

Mainline VIDE contributions from Ben Wandelt, Nico Hamaus, Alice Pisani, 
Paul Zivick, and Qingqing Mao.
This toolkit includes ZOBOV, originally developed by Mark Neyrinck. 
See `zobov/zobov_readme.txt` for copyright/license information. 
SDF library provided by Michael S. Warren and John Salmon. 
HOD fitting code provided by Francisco Navarro. 
HOD halo population code provided by Jeremy Tinker.
RAMSES module provided by Benjamin B. Thompson.


Requirements
------------

The package swig needs to be installed and available in the PATH (http://www.swig.org/). It 
is required by scipy and we have not decided to bundle it with VIDE at the moment.


Quick Start Guide
-----------------

After installing the package with `python3 setup.py install --user`, you can execute

```
python3 -m void_pipeline  your_config_file.py
```

The VIDE tools are all packaged in the `vide` package. 



Notes for CONDA
---------------


If you use a conda installation, you have to be sure to use all the building tools that 
are consistent. On linux that means for example installing the conda packages `gcc_linux-64`
 and `gxx_linux-64`. In addition to that it is recommended to define the environment variable 
`LIBRARY_PATH=the_path_of_your_conda_environment_with_/lib`. For example if your environment
is in '/home/user/conda' you should define

```
export LIBRARY_PATH=/home/user/conda/lib
```

You can then initiate the construction with

```
python3 setup.py build
```

Version Summary
-----------------

v1.0 - Initial Release
v2.0 - Ported to python3, revisited build system
