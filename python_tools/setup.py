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
       ['void_python_tools','void_python_tools.backend','void_python_tools.apTools',
        'void_python_tools.apTools.profiles','void_python_tools.apTools.chi2', 'void_python_tools.plotting'],
    #ext_modules = [Extension("void_python_tools.chi2.velocityProfileFitNative", ["void_python_tools/chi2/velocityProfileFitNative.pyx"], libraries=["gsl", "gslcblas"]), Extension("void_python_tools.chi2.likelihoo", ["void_python_tools/chi2/likelihood.pyx"], libraries=["gsl", "gslcblas"])]
    ext_modules = [
       Extension("void_python_tools.apTools.chi2.velocityProfileFitNative", 
                 ["void_python_tools/apTools/chi2/velocityProfileFitNative.pyx"],
       libraries=["gsl", "gslcblas"], library_dirs=[VOID_GSL+"/lib"], include_dirs=[VOID_GSL+"/include"])
    ]
)
