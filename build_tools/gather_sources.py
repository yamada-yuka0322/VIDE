#+
#   VIDE -- Void IDentification and Examination -- ./build_tools/gather_sources.py
#   Copyright (C) 2010-2014 Guilhem Lavaux
#   Copyright (C) 2011-2014 P. M. Sutter
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; version 2 of the License.
# 
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program; if not, write to the Free Software Foundation, Inc.,
#   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#+
import shutil
import tempfile
import re
from git import Repo,Tree,Blob

def apply_license(license, relimit, filename):
  header = re.sub(r'@FILENAME@', filename, license)

  with open(filename) as f:
    lines = f.read()

  lines = re.sub(relimit, lambda x: (("" if x.group(1) is None else x.group(1)) + header), lines)

  with tempfile.NamedTemporaryFile(delete=False, mode="wt", encoding="UTF-8") as tmp_sources:
    tmp_sources.write(lines)

  shutil.move(tmp_sources.name, filename)


def apply_python_license(filename):
  license="""#+
#   VIDE -- Void IDentification and Examination -- @FILENAME@
#   Copyright (C) 2010-2014 Guilhem Lavaux
#   Copyright (C) 2011-2014 P. M. Sutter
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; version 2 of the License.
# 
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program; if not, write to the Free Software Foundation, Inc.,
#   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#+
"""
   
  print("Shell/Python file: %s" % filename)
  relimit=r'^(#!.*\n)?#\+\n(#.*\n)*#\+\n'
  apply_license(license, relimit, filename)


def apply_cpp_license(filename):
  license="""/*+
    VIDE -- Void IDentification and Examination -- @FILENAME@
    Copyright (C) 2010-2014 Guilhem Lavaux
    Copyright (C) 2011-2014 P. M. Sutter

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; version 2 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
+*/
"""
  print("C++ file: %s" % filename)
  relimit = r'(?s)^()/\*\+.*\+\*/\n'
  apply_license(license, relimit, filename)
  

def analyze_tree(prefix, t):
  for entry in t:
    ename = entry.name
    if ename == 'external' or ename == 'zobov':
      continue
    if type(entry) == Tree:
      analyze_tree(prefix + "/" + ename, entry)
    elif type(entry) == Blob:
      if re.match(".*\.(sh|py|pyx)$",ename) != None:
        fname=prefix+"/"+ename
        apply_python_license(fname)
      if re.match('.*\.(cpp|hpp|h)$', ename) != None:
        fname=prefix+"/"+ename
        apply_cpp_license(fname)


if __name__=="__main__":
  import sys
  if len(sys.argv) >= 2:
      apply_python_license(sys.argv[1])
      sys.exit(0)

  repo = Repo(".")
  assert repo.bare == False
  t = repo.tree()
  analyze_tree(".", t)
