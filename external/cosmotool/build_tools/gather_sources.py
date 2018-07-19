#+
# This is CosmoTool (./build_tools/gather_sources.py) -- Copyright (C) Guilhem Lavaux (2007-2014)
#
# guilhem.lavaux@gmail.com
#
# This software is a computer program whose purpose is to provide a toolbox for cosmological
# data analysis (e.g. filters, generalized Fourier transforms, power spectra, ...)
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#+
import shutil
import tempfile
import re
from git import Repo,Tree,Blob

def apply_license(license, relimit, filename):
  header = re.sub(r'@FILENAME@', filename, license)

  f = file(filename)
  lines = f.read()
  f.close()

  lines = re.sub(relimit, '', lines)
  lines = header + lines

  with tempfile.NamedTemporaryFile(delete=False) as tmp_sources:
    tmp_sources.write(lines)

  shutil.move(tmp_sources.name, filename)


def apply_python_license(filename):
  license="""#+
# This is CosmoTool (@FILENAME@) -- Copyright (C) Guilhem Lavaux (2007-2014)
#
# guilhem.lavaux@gmail.com
#
# This software is a computer program whose purpose is to provide a toolbox for cosmological
# data analysis (e.g. filters, generalized Fourier transforms, power spectra, ...)
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#+
"""
   
  print("Shell/Python file: %s" % filename)
  relimit=r'^#\+\n(#.*\n)*#\+\n'
  apply_license(license, relimit, filename)


def apply_cpp_license(filename):
  license="""/*+
This is CosmoTool (@FILENAME@) -- Copyright (C) Guilhem Lavaux (2007-2014)

guilhem.lavaux@gmail.com

This software is a computer program whose purpose is to provide a toolbox for cosmological
data analysis (e.g. filters, generalized Fourier transforms, power spectra, ...)

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
+*/
"""
  relimit = r'^(?s)/\*\+.*\+\*/\n'
  print("C++ file: %s" % filename)
  apply_license(license, relimit, filename)
  

def analyze_tree(prefix, t):
  for entry in t:
    ename = entry.name
    if ename == 'external' or ename == 'zobov':
      continue
    if type(entry) == Submodule:
      continue
    elif type(entry) == Tree:
      analyze_tree(prefix + "/" + ename, entry)
    elif type(entry) == Blob:
      if ename == './src/hdf5_flash.h' or ename == './src/h5_readFlash.cpp' or ename == './src/h5_readFlash.hpp':
        continue

      if re.match(".*\.(sh|py|pyx)$",ename) != None:
        fname=prefix+"/"+ename
        apply_python_license(fname)
      if re.match('.*\.(cpp|hpp|h)$', ename) != None:
        fname=prefix+"/"+ename
        apply_cpp_license(fname)


if __name__=="__main__":
  repo = Repo(".")
  assert repo.bare == False
  t = repo.tree()
  analyze_tree(".", t)
