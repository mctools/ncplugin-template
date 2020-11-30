#!/bin/bash

set -e
set -u

if [ "x${1:-dummy}" == x-h -o "x${1:-dummy}" == x--help ]; then
    echo 'Usage:'
    echo
    echo "1) Run with no arguments to download, build, and install latest commit in NCrystal master branch on github:"
    echo
    echo "  ${0}"
    echo
    echo "2) Download, build, and install specific commit, branch, or tag from NCrystal repo at github. For example:"
    echo
    echo "  ${0} gitref develop     #<--- latest commit in develop branch"
    echo "  ${0} gitref v2.2.0      #<--- specific NCrystal release via tag"
    echo
    echo "1) Unpack, build, and install from tar-ball with ncrystal sources:"
    echo
    echo "  ${0} some-downloaded-ncrystal-sources.tar.gz"
    echo
    exit 0
fi

export _ncputildir="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [ ! -f "${_ncputildir}/../../ncplugin_name.txt" ]; then
    echo "ERROR: Script not in expected location..."
    exit 1
fi

export ncbase="${_ncputildir}/cache/ncrystal"
gitref=master

if [[ "${1:-dummy}" == gitref ]]; then
    gitref="${2}"
    shift 2
fi


if [ -e "${ncbase}" ]; then
    echo "ERROR: Directory ${ncbase} already exists. Please remove it and rerun."
    exit 1
fi

for cmd in curl tar cmake python3; do
    which $cmd >/dev/null 2>&1 || ( echo "ERROR: Required command not found: $cmd" ; exit 1 )
done

python3 "${_ncputildir}/validateprereqs.py" --skip-ncrystal-check

#Build (respect CMAKE_BUILD_PARALLEL_LEVEL if set, otherwise detect):
export CMAKE_BUILD_PARALLEL_LEVEL=${CMAKE_BUILD_PARALLEL_LEVEL:-$(python3 -c 'import os;print(os.cpu_count())')}
if [ $? != 0 ]; then
    echo "ERROR could not determine number of processes to use for build. Please set CMAKE_BUILD_PARALLEL_LEVEL env var and rerun"
    return 1
fi

mkdir -p "${ncbase}"
cd ${ncbase}

dlcmd="curl -LsS https://api.github.com/repos/mctools/ncrystal/tarball/${gitref}  --output src.tar.gz"
srcdirpattern="*ncrystal*"
if [[ "${1:-dummy}" == *.tar.gz ]]; then
    dlcmd="cp ${1} src.tar.gz"
    gitref='n/a'
    shift 1
fi
echo "Attempting to download latest NCrystal release (git ref: ${gitref})"
echo "Download command is: $dlcmd"
mkdir download
cd download
$dlcmd
cd ../
tar xf download/src.tar.gz
mv ${srcdirpattern}/ src/
rm -f src.tar.gz
echo "Downloaded and unpacked NCrystal v$(cat src/VERSION)"
echo "Source dir is: ${ncbase}/src"
mkdir bld
mkdir install
baseflags="-DBUILD_G4HOOKS=OFF -DBUILD_EXAMPLES=OFF -DBUILD_EXTRA=OFF -DCMAKE_BUILD_TYPE=RelWithDebInfo"
echo "Configuring NCrystal."
echo "   Applying flags: ${baseflags}"
echo "   Applying extra flags passed on command line:" "$@"
cmake -S "${ncbase}/src" -B "${ncbase}/bld" -DCMAKE_INSTALL_PREFIX="${ncbase}/install" ${baseflags} "$@" \
    || ( echo "Aborting due to configuration errors"; exit 1 )
echo "Building NCrystal with $CMAKE_BUILD_PARALLEL_LEVEL processes."
cmake --build "${ncbase}/bld" || ( echo "Aborting due to build errors."; exit 1 )
echo "Installing NCrystal to ${ncbase}/install"
cmake --install "${ncbase}/bld" || ( echo "Aborting due to errors during installation."; exit 1 )
echo
echo "Done. Installed NCrystal v$(cat src/VERSION) (git ref ${gitref})."
echo "Location: ${ncbase}/install"
echo
echo "If you wish to save space, you can remove the build leftovers with:"
echo "   rm -rf ${ncbase}/bld"
echo
echo "This NCrystal installation will be automatically activated when you source testcode/util/bootstrap.sh"
