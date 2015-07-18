#BEAR: Balance Equations for Atomic Reactions (in dev)

## License 
BEAR  is distributed under the terms of the GNU Lesser General Public Licence version 3 (LGPL) version 3.

##Getting started


### Step by Step installation

#### 1. Build Requirements

cmake >= 2.8.2

C++11 Compiler:

* GCC >= 4.6
* clang >= 3.2

Third party libraries:
* boost>=1.41 (required)
* boost-numeric-bindings and Lapack>=3.5 (required for solving differential systems)
* Root>=5.32 (optional)

#### 2. Install [BEAR](https://github.com/NicolasWinckler/BEAR)

    "Bear_path_install" is used here as a directory name, but can be renamed as you want.
    
    ```bash
    # Set the shell variable SIMPATH to the boost installation directory
    export SIMPATH=~/boost_install
    [setenv SIMPATH ~/boost_install]

    cd ~/Bear_path_install
    git clone https://github.com/NicolasWinckler/BEAR
    cd BEAR
    mkdir build
    cd build
    make
    make install
    ```
