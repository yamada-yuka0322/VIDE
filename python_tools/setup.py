#+
#   VIDE -- Void IDentification and Examination -- ./python_tools/setup.py
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
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np
import os

VOID_GSL=os.environ.get('VOID_GSL')

setup(
    name='void_python_tools',
    version='1.0',
    cmdclass = {'build_ext': build_ext},
    include_dirs = [np.get_include()],
    packages=
       ['void_python_tools','void_python_tools.backend','void_python_tools.apTools', 'void_python_tools.xcor', 'void_python_tools.voidUtil',
        'void_python_tools.apTools.profiles','void_python_tools.apTools.chi2',],
    #ext_modules = [Extension("void_python_tools.chi2.velocityProfileFitNative", ["void_python_tools/chi2/velocityProfileFitNative.pyx"], libraries=["gsl", "gslcblas"]), Extension("void_python_tools.chi2.likelihoo", ["void_python_tools/chi2/likelihood.pyx"], libraries=["gsl", "gslcblas"])]
    #ext_modules = [
    #   Extension("void_python_tools.apTools.chi2.velocityProfileFitNative", 
    #             ["void_python_tools/apTools/chi2/velocityProfileFitNative.pyx"],
    #   libraries=["gsl", "gslcblas"], library_dirs=[VOID_GSL+"/lib"], include_dirs=[VOID_GSL+"/include"])
    #]
)
