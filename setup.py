#+
#   VIDE -- Void IDentification and Examination -- ./setup.py
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
import stat
import os
import sys
import shutil
from distutils.command.install_data import install_data
from distutils.command.build_scripts import build_scripts
import pathlib
from setuptools import find_packages, setup, Extension, Command
from setuptools.command.build_ext import build_ext
from setuptools.command.install_lib import install_lib
from setuptools.command.install_scripts import install_scripts
import struct

BITS = struct.calcsize("P") * 8
PACKAGE_NAME = "vide"

class CMakeExtension(Extension):
    """
    An extension to run the cmake build

    This simply overrides the base extension class so that setuptools
    doesn't try to build your sources for you
    """

    def __init__(self, name, sources=["dummy_extension/empty.c"]):

        super().__init__(name = name, sources = sources)

        self.SOURCE_DIR = str(pathlib.Path().absolute())

class InstallCMakeLibsData(install_data):
    """
    Just a wrapper to get the install data into the egg-info

    Listing the installed files in the egg-info guarantees that
    all of the package files will be uninstalled when the user
    uninstalls your package through pip
    """

    def run(self):
        """
        Outfiles are the libraries that were built using cmake
        """

        # There seems to be no other way to do this; I tried listing the
        # libraries during the execution of the InstallCMakeLibs.run() but
        # setuptools never tracked them, seems like setuptools wants to
        # track the libraries through package data more than anything...
        # help would be appriciated

        self.outfiles = self.distribution.data_files

class InstallCMakeLibs(install_lib):
    """
    Get the libraries from the parent distribution, use those as the outfiles

    Skip building anything; everything is already built, forward libraries to
    the installation step
    """

    def run(self):
        """
        Copy libraries from the bin directory and place them as appropriate
        """

        self.announce("Moving library files", level=3)

        self.distribution.run_command("install_data")
        self.distribution.run_command("install_scripts")

        super().run()

class BuildCMakeExt(build_ext):
    """
    Builds using cmake instead of the python setuptools implicit build
    """

    def run(self):
        """
        Perform build_cmake before doing the 'normal' stuff
        """

        for extension in self.extensions:

            if extension.name == 'vide':
                self.package = 'vide'

                self.build_cmake(extension)

        super().run()

    def build_cmake(self, extension: Extension):
        """
        The steps required to build the extension
        """

        self.announce("Preparing the build environment", level=3)

        package_dir = os.path.abspath(os.path.join(self.build_lib, 'vide'))


        extension.build_dir = pathlib.Path(self.build_temp)
        extension.bin_dir = str(pathlib.Path(os.path.join(extension.build_dir, 'private_install')).absolute())
        SOURCE_DIR = extension.SOURCE_DIR
        build_dir = extension.build_dir

        extension_path = pathlib.Path(self.get_ext_fullpath(extension.name))

        os.makedirs(build_dir, exist_ok=True)
        os.makedirs(extension_path.parent.absolute(), exist_ok=True)

        cython_code = os.path.join(str(build_dir.absolute()),'mycython')
        with open(cython_code, mode="wt") as ff:
          ff.write(f"#!{sys.executable}\n"
                    "from Cython.Compiler.Main import setuptools_main\n"
                    "setuptools_main()")
        os.chmod(cython_code, stat.S_IXUSR|stat.S_IWUSR|stat.S_IRUSR|stat.S_IRGRP)

        # Now that the necessary directories are created, build

        self.announce("Configuring cmake project", level=3)

        # Change your cmake arguments below as necessary
        # Below is just an example set of arguments for building Blender as a Python module

        PYTHON_bin_package=f"{build_dir.absolute()}/private_install"

        c_compiler=os.environ.get('CC', get_config_var("CC"))
        cxx_compiler=os.environ.get('CXX', get_config_var("CXX"))

        self.spawn(['cmake', '-H'+SOURCE_DIR, '-B'+self.build_temp,
                    '-DENABLE_OPENMP=ON','-DINTERNAL_BOOST=ON','-DINTERNAL_EIGEN=ON',
                    '-DINTERNAL_HDF5=ON','-DINTERNAL_NETCDF=ON','-DINTERNAL_GSL=ON',
                    '-DBUILD_PYTHON=ON', '-DINSTALL_PYTHON_LOCAL=OFF',
                    '-DCOSMOTOOL_PYTHON_PACKAGING=ON',
                    f"-DCYTHON={cython_code}",
                    '-DINSTALL_CTOOLS_IN_PYTHON=ON',
                    f"-DCMAKE_C_COMPILER={c_compiler}", f"-DCMAKE_CXX_COMPILER={cxx_compiler}",
                    f"-DPYTHON_SITE_PACKAGES={PYTHON_bin_package}",
                    f"-DPYTHON_EXECUTABLE={sys.executable}"])


        self.announce("Building binaries", level=3)

        self.spawn(["cmake", "--build", self.build_temp, "--target", "all",
                      "--config", "Release","--","VERBOSE=1"])

        self.spawn(["cmake", "--build", self.build_temp, "--target", "install",
                      "--config", "Release","--","VERBOSE=1"])

        # Build finished, now copy the files into the copy directory
        # The copy directory is the parent directory of the extension (.pyd)

        self.announce("Moving built python module", level=3)

        bin_dir = PYTHON_bin_package
        #extension.bin_dir
        self.distribution.bin_dir = bin_dir
        target_dir = os.path.abspath(os.path.join(package_dir,'bin'))
        print(target_dir)
        shutil.rmtree(target_dir, ignore_errors=True)
        shutil.move(f"{PYTHON_bin_package}/void_python_tools/bin", target_dir)

        shutil.copy(f"{self.build_temp}/pipeline/prepareInputs.py", f"{self.build_temp}/vide_prepare_simulation")

class VideScripts(build_scripts):

    user_options = build_scripts.user_options + [
        ('build-temp=', 't', "temporary directory where scripts are stored"),
    ]

    vide_scripts = ["vide_prepare_simulation"]

    def initialize_options(self):
        super(VideScripts, self).initialize_options()
        self.build_temp = None

    def finalize_options(self):
        super(VideScripts, self).finalize_options()
        self.set_undefined_options('build_ext',
                                   ('build_temp', 'build_temp'))
        self.scripts = [os.path.join(self.build_temp, v) for v in self.vide_scripts]

    def run(self):
        self.copy_scripts()

vide_extension = CMakeExtension(name="vide")

setup(name='vide',
      version='2.0',
      packages=find_packages('python_tools'),
      package_dir={'': 'python_tools'},
      setup_requires=['cython','setuptools','healpy','argparse','scipy','astropy','extension-helpers','netCDF4'],
      install_requires=['argparse','scipy','astropy','extension-helpers','netCDF4','healpy'],
      ext_modules=[vide_extension],
      description='The VIDE pipeline analysis for Cosmic Voids',
      long_description=open("./README.md", 'r').read(),
      long_description_content_type="text/markdown",
      keywords="cosmology, interpolation, cmake, extension",
      classifiers=["Intended Audience :: Developers",
                   "License :: OSI Approved :: "
                   "GNU Lesser General Public License v3 (LGPLv3)",
                   "Natural Language :: English",
                   "Programming Language :: C",
                   "Programming Language :: C++",
                   "Programming Language :: Python",
                   "Programming Language :: Python :: 3",
                   "Programming Language :: Python :: Implementation :: CPython"],
      license='CeCILL-v2',
      scripts=["vide_prepare_simulation"],
      cmdclass={
          'build_ext': BuildCMakeExt,
          'install_data': InstallCMakeLibsData,
          'install_lib': InstallCMakeLibs,
          'build_scripts': VideScripts,
          'install_scripts': install_scripts,
          }
)
