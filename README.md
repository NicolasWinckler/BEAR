#BEAR: Balance Equations for Atomic Reactions 

## License 
BEAR  is distributed under the terms of the GNU Lesser General Public Licence version 3 (LGPL) version 3.

##Getting started


### Step by Step installation

#### 1. Build Requirements

cmake >= 2.8.3

C++11 Compiler:

* GCC >= 4.6
* clang >= 3.2

Third party libraries:
* Boost>=1.41 (required)
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
