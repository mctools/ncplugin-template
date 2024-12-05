# NCrystal plugin CrysExtn

This plugin is dedicated to study the extinction effects observed in polycrystalline matherials. Both Sabine's and Becker and Coppens' (BC) extinction models are implemented.

# Installation

The installation necessitates to recompile NCrystal, a tutorial is presented below:

```
conda create -n ncrystal-extinction -y -c conda-forge jupyter numpy scipy matplotlib pip h5py
conda activate ncrystal-extinction
cd path/to/directory
git clone https://github.com/mctools/ncrystal.git
git clone https://github.com/XuShuqi7/ncplugin-CrysExtn.git
cd ncrystal
mkdir build
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DNCRYSTAL_BUILTIN_PLUGINS="path/to/directory/ncplugin-CrysExtn"
make -j 8
make install
$($CONDA_PREFIX/bin/ncrystal-config --setup)
nctool --test
python3 -c 'import NCrystal; NCrystal.test()'
python3 -c 'import NCrystal; print(NCrystal.__version__)'
```

# Note

A jupyter notebook is attached in the folder "verification" for testing the implementation. The extinction model will be further implemented inside NCrystal as an internal funtionality, please check the lastest version of the code: https://github.com/mctools/ncrystal
