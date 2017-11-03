# How to build OpenMM from source

This recipe build OpenMM from source and install its Python API into a separated conda environment.

```bash
# FIRST MAKE SURE THE OLD BUILD IS UNINSTALLED (SEE NEXT SECION).

# Prepare and install dependencies (everything else should be included in OS X)
conda create -n testopenmm swig numpy  # numpy is needed for the make PythonInstall step
source activate testopenmm
conda install -c dlr-sc doxygen  # there is also an omnia doxygen package but not for osx64
git clone git@github.com:pandegroup/openmm.git openmm.pandegroup

# If you want to have easy life and browse/edit the source code with Xcode
mkdir openmm_xcodeproj
cd openmm_xcodeproj
cmake -G Xcode ../openmm.pandegroup/  # automatically creates XCode project file
cd ..

# Build and install source code
mkdir build_openmm
cd build_openmm
# See http://docs.openmm.org/latest/userguide/library.html#all-platforms for options.
# I kept /usr/local/openmm/ as installation path but you can change it if you don't have root privileges.
# Make sure the python interpreter is set to the one in miniconda/envs/testopenmm/bin/python.
ccmake -i ../openmm.pandegroup/
make
sudo make install

# Install the Python API in your testopenmm conda environment
source activate testopenmm  # if you haven't done it before
make PythonInstall
```

## To uninstall
```bash
conda remove --name testopenmm --all
conda clean -itps --yes
sudo rm -r /usr/local/openmm/
```

## Tested on
Issues found when compiling C++ code can vary widely among different systems/compilers. This recipe worked on (add yours)
* Macbook pro, OS X 10.10.5, clang 7.0.2

