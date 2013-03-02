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
#   VIDE -- Void IDEntification pipeline -- @FILENAME@
#   Copyright (C) 2010-2013 Guilhem Lavaux
#   Copyright (C) 2011-2013 Paul M. Sutter
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
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
  relimit=r'^#\+\n(#.*\n)*#\+\n'
  apply_license(license, relimit, filename)


def apply_cpp_license(filename):
  license="""/*+
    VIDE -- Void IDEntification pipeline -- @FILENAME@
    Copyright (C) 2010-2013 Guilhem Lavaux
    Copyright (C) 2011-2013 Paul M. Sutter

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
+*/
"""
  relimit = r'^(?s)/\*\+.*\+\*/'
  print("C++ file: %s" % filename)
  apply_license(license, relimit, filename)
  

def analyze_tree(prefix, t):
  for ename,entry in t.items():
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
  repo = Repo(".")
  assert repo.bare == False
  t = repo.tree()
  analyze_tree(".", t)
