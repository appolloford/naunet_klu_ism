# naunet_klu_ism

This is a C++ implementation of the chemical network used in [Walsh et al. 2015](https://ui.adsabs.harvard.edu/abs/2015A%26A...582A..88W/abstract). 
It is also a test example of [naunet](https://github.com/appolloford/naunet).

## Dependencies

- [suitesparse](https://github.com/DrTimothyAldenDavis/SuiteSparse)^5.7.2 (tested by sundials)
- [sundials](https://github.com/LLNL/sundials)^5.6.1

### Install denpendencies

As the KLU solver in sundials requires the suitesparse library, suitesparese must be installed before sundials.

- MacOS: 
  1. `brew install suite-sparse`
  2. `brew` supported installing sundials, but it was updated to 6.2.0 recently. Here is an example showing how to install from source
     ```bash
     $ wget https://github.com/LLNL/sundials/releases/download/v5.7.0/sundials-5.7.0.tar.gz
     $ tar -zxvf sundials-5.7.0.tar.gz
     $ mkdir build-sundials && cd build-sundials
     $ cmake ../sundials-5.7.0 -DCMAKE_INSTALL_PREFIX=/usr/local/sundials \
              -DSUNDIALS_INDEX_SIZE=32 -DENABLE_KLU=ON \
              -DKLU_INCLUDE_DIR=$(brew --prefix suite-sparse)/include \
              -DKLU_LIBRARY_DIR=$(brew --prefix suite-sparse)/lib \
              -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_INSTALL=OFF
     $ make -j4 && sudo make install
     ```

- Ubuntu:
  1. `suitesparse` can be install by `apt-get`, I did not test it because the version is earlier, but you can try.
     Here is an example of installing from sources. It may take for a while, or you can install only the required part.
     ```bash
     $ sudo apt-get install libgmp3-dev libmpc-dev # install suitesparse dependencies
     $ wget https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v5.9.0.tar.gz 
     $ tar -zxvf v5.9.0.tar.gz && cd SuiteSparse-5.9.0
     $ make library -j4 && sudo make install
     ```

## Build & Test

1. Clone this repository: `git clone https://github.com/appolloford/naunet_klu_ism.git`
2. Go into the folder: `cd naunet_klu_ism`
3. Configure the project: `cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release`. 
   You may need to add `-DSUNDIALS_DIR=/path/to/sundials/lib/cmake/sundials` if sundials is not installed in the system path.
   If you want to use python module add `-DMAKE_PYTHON` in the end.
   To install the project as a library, add `-DCMAKE_INSTALL_PREFIX=.` in the end.
4. Build project: `cmake --build build`
5. To test, go into `build` folder and run `ctest`: `cd build && ctest`. The executables are in `build/tests/`
6. If you want to install it, execute `make install` in `build` folder

## Usage

TBD. (Check the examples in `tests/`)
