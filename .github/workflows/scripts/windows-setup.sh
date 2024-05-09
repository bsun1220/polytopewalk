#!/bin/bash

# Clone and bootstrap vcpkg
git clone https://github.com/microsoft/vcpkg
./vcpkg/bootstrap-vcpkg.bat
./vcpkg/vcpkg integrate install

# Install eigen3 and ipopt with vcpkg
./vcpkg/vcpkg install eigen3
./vcpkg/vcpkg install coin-or-ipopt

# Install ifopt from source
git clone https://github.com/ethz-adrl/ifopt.git && cd ifopt
mkdir build && cd build
cmake ..
make
make install
