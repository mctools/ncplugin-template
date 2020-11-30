#!/usr/bin/env python3
import sys
if sys.version_info[0:2] < (3,6):
    raise SystemExit('Error: python 3.6 or later required')


import pathlib
import os
import shutil
import subprocess

cmake_min_version=(3,15)
if not shutil.which('cmake'):
    raise SystemExit("Error: cmake command not found. Please install CMake (at least version %i.%i)."%cmake_min_version)
cmake_version_str = subprocess.check_output(["cmake", "--version"]).decode().split()[2]
cmake_version = tuple(int(e) for e in cmake_version_str.split(".")[0:2])
if cmake_version < cmake_min_version:
    raise SystemExit("CMake version is too low. Please use at least CMake %i.%i or later."%cmake_min_version)

if not '--skip-ncrystal-check' in sys.argv[1:]:
    if not shutil.which('ncrystal-config'):
        raise SystemExit("Error: ncrystal-config command not found. Do you have an active and updated NCrystal installation present? If not you can try to install by running the command ./testcode/utils/installncrystal.x")

for f in pathlib.Path(__file__).parent.glob('../scripts/*'):
    if any(c in f.name for c in '~#. '):
        continue#ignore backup files etc.
    if not os.access(f, os.X_OK):
        raise SystemExit(f'Script file {f.name} is not executable. Please rectify by running the command:\n\n  chmod +x {f.absolute().resolve()}\n')
