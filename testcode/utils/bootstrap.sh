#!/bin/bash

export _ncputildir="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [ ! -f "${_ncputildir}/../../ncplugin_name.txt" ]; then
    echo "ERROR: Script not in expected location..."
    return 1
fi

which ncrystal-config 2>&1 >/dev/null
if [ $? != 0 -a -f "${_ncputildir}/cache/ncrystal/install/setup.sh" ]; then
    echo "Activating existing ncrystal installation at ${_ncputildir}/cache/ncrystal/install"
    . "${_ncputildir}/cache/ncrystal/install/setup.sh"
fi

python3 "${_ncputildir}"/validateprereqs.py
if [ $? != 0 ]; then
    echo
    echo "ERROR: Basic prerequisites check failed."
    echo
    return 1
fi

export _ncpsrcdir="$( cd -P "${_ncputildir}/../../" && pwd )"
export _ncpblddir="${_ncputildir}/cache/bld"
export _ncpinstdir="${_ncputildir}/cache/install"

function ncpluginbuild() {
    #Check options:
    local _ncpdoforce
    local _ncpdoclean
    _ncpdoforce=0
    _ncpdoclean=0
    if [ x"${1:-}" == x--force -o x"${1:-}" == x-f ]; then
        _ncpdoforce=1
        shift 1
    elif [ x"${1:-}" == x--clean -o x"${1:-}" == x-c ]; then
        _ncpdoclean=1
        shift 1
    elif [ x"${1:-}" == x--help -o x"${1:-}" == x-h ]; then
        local usage
        usage=$(cat <<-END
The ncpluginbuild command takes care of calling CMake to configure / build /
install the plugin and testcode, and subsequently set up environment variables
so the scripts and C++ applications can be invoked from the command line, the
Python module imported, and the libraries loaded, as needed.

Typical usage is simply to invoke the ncpluginbuild command with no arguments,
and to invoke it again whenever some C++ code was changed or a python module or
script file was added or removed. More detailed usage instructions are available
below.

Note that builds will be default try to use all available processes on your
machine.  To override this, you can set the CMAKE_BUILD_PARALLEL_LEVEL
environment variable to the desired number of processes to use.

Usage:

  ncpluginbuild [--force|-f] <CMAKEARGS>

                             : Configure, build and install both plugin and
                               testcode, and modify environment so all resulting
                               scripts and C++ applications are available to run
                               afterwards. Note that python modules and scripts
                               are symlinked into the installation area, so if you
                               are just changing such files, you do not need to
                               rerun the ncpluginbuild command in order for the
                               changes to take effect.

                               If any <CMAKEARGS> are provided, they will be
                               passed along to CMake for the configuration
                               step. If running again with different
                               <CMAKEARGS>, a full reconfiguration and rebuild
                               will automatically be performed from scratch.

                               The --force or -f flag can be used to force a
                               full reconfiguration and rebuild (as opposed to
                               just an incremental build). Ideally, this should
                               never be needed, but it is here as a fall-back
                               option in case something does not work as expected.

  ncpluginbuild (--clear|-c) : clean bld output, environment variables and even
                               delete the ncpluginbuild command. (you must
                               ". testcode/util/bootstrap.sh" to get it back
                               again).

  ncpluginbuild (--help|-h)  : Print these usage instructions
END
)
        echo
        echo "${usage}"
        echo
        return
    fi
    python3 "${_ncputildir}"/validateprereqs.py
    if [ $? != 0 ]; then
        echo "ncpluginbuild: Aborting due to errors."
        return 1
    fi
    #Clear installation area every time (also clear any references to this
    #directory in the environment):
    echo "ncpluginbuild: clearing ${_ncpinstdir}"
    rm -rf "${_ncpinstdir}"
    function ncplugin_prunepath() {
        P=$(IFS=:;for p in ${!1:-}; do [[ $p != ${2}* ]] && echo -n ":$p" || :; done)
        export $1="${P:1:99999}"
    }
    ncplugin_prunepath PATH "${_ncpinstdir}"
    ncplugin_prunepath LD_LIBRARY_PATH "${_ncpinstdir}"
    ncplugin_prunepath DYLD_LIBRARY_PATH "${_ncpinstdir}"
    ncplugin_prunepath PYTHONPATH "${_ncpinstdir}"
    ncplugin_prunepath NCRYSTAL_DATA_PATH "${_ncpinstdir}"
    unset ncplugin_prunepath
    if [ $_ncpdoclean == 1 ]; then
        rm -rf "${_ncpblddir}"
        unset _ncputildir _ncpsrcdir _ncpblddir _ncpinstdir
        unset ncpluginbuild
        echo "ncpluginbuild: Cleared environment variables, removed bld directories and the ncpluginbuild command."
        echo "ncpluginbuild: Your must source the bootstrap.sh file again when/if you wish to resume working with the ncpluginbuild command."
        return 0
    fi
    if [ $_ncpdoforce == 0 -a -d "${_ncpblddir}" -a ! -f "${_ncpblddir}"/CMakeCache.txt ]; then
        echo "ncpluginbuild: Missing cache. Auto-enabling --force mode."
        _ncpdoforce=1
    fi
    local oldargs
    local newargs
    oldargs=$(cat ${_ncpblddir}/ncplugin_usercmakeargs.txt 2>/dev/null || echo "..<n/a>..")
    newargs="$(echo "$@")"
    if [ -f "${_ncpblddir}"/CMakeCache.txt -a x"$oldargs" != x"$newargs" ]; then
        echo "ncpluginbuild: Argument changes detected. Auto-enabling --force mode."
        _ncpdoforce=1
    fi
    if [ $_ncpdoforce == 1 ]; then
        echo "ncpluginbuild: (--force mode) clearing ${_ncpblddir}"
        rm -rf "${_ncpblddir}"
    fi
    if [ ! -d "${_ncpblddir}" ]; then
        echo "ncpluginbuild: configuring with build dir ${_ncpblddir}"
        if [ $# != 0 ]; then
            echo "ncpluginbuild: Using extra CMake configuration arguments:" "$@"
        fi
        mkdir -p "${_ncpblddir}" && \
            cmake -S "${_ncpsrcdir}" -B "${_ncpblddir}" \
                  -DCMAKE_INSTALL_PREFIX="${_ncpinstdir}" -DNCPLUGIN_INSTALLSYMLINKS=ON "$@"
        if [ $? != 0 ]; then
            echo "ncpluginbuild: Aborting due to configuration errors"
            return 1
        fi
        echo -n "$newargs" > ${_ncpblddir}/ncplugin_usercmakeargs.txt
    fi
    (
        #Build (respect CMAKE_BUILD_PARALLEL_LEVEL if set, otherwise detect):
        export CMAKE_BUILD_PARALLEL_LEVEL=${CMAKE_BUILD_PARALLEL_LEVEL:-$(python3 -c "import os;print(os.cpu_count())")}
        if [ $? != 0 ]; then
            echo "ncpluginbuild: ERROR could not determine number of processes to use for build. Please set CMAKE_BUILD_PARALLEL_LEVEL env var and rerun"
            return 1
        fi
        echo "ncpluginbuild: Building with $CMAKE_BUILD_PARALLEL_LEVEL processes."
        cmake --build "${_ncpblddir}"
        if [ $? != 0 ]; then
            echo "ncpluginbuild: Aborting due to build errors (if weird try to rerun with --force flag)."
            return 1
        fi
    ) || return 1
    (
        echo "ncpluginbuild: installing to  ${_ncpinstdir}"
        cmake --install "${_ncpblddir}"
        if [ $? != 0 ]; then
            echo "ncpluginbuild: Aborting due to errors during installation."
            return 1
        fi
    ) || return 1
    echo "ncpluginbuild: Installation completed, updating environment variables to activate."
    export PATH="${_ncpinstdir}/bin:${PATH:-}"
    export LD_LIBRARY_PATH="${_ncpinstdir}/lib:${LD_LIBRARY_PATH:-}"
    export DYLD_LIBRARY_PATH="${_ncpinstdir}/lib:${DYLD_LIBRARY_PATH:-}"
    export PYTHONPATH="$(echo "${_ncpinstdir}"/share/NCPluginTestCode_*/python):${PYTHONPATH:-}"
    export NCRYSTAL_DATA_PATH="$(echo "${_ncpinstdir}"/share/NCPluginTestCode_*/data):${NCRYSTAL_DATA_PATH:-}"
    echo
    echo "ncpluginbuild: Done! You can run your scripts/applications by invoking the relevant command from the list below:"
    echo
    (
        cd "${_ncpinstdir}/bin"
        for cmd in *; do
            echo "    $> $cmd"
        done
    )
    echo
}

echo "Done. You now have the ncpluginbuild command available. Run \"ncpluginbuild --help\" to get started."
