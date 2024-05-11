#!/bin/bash

# Clone and bootstrap vcpkg
git clone https://github.com/microsoft/vcpkg
./vcpkg/bootstrap-vcpkg.bat
./vcpkg/vcpkg integrate install

# Install eigen3 and ipopt with vcpkg
./vcpkg/vcpkg install eigen3
./vcpkg/vcpkg install pkgconf
# ./vcpkg/vcpkg install coin-or-ipopt

# Test ipopt install from source
# Install Mumps
git clone https://github.com/coin-or-tools/ThirdParty-Mumps.git && cd ThirdParty-Mumps
./get.Mumps
./configure
make; make install
cd ..

git clone https://github.com/coin-or/Ipopt.git && cd Ipopt
./configure
make; make test; make install
cd ..

# Install ifopt from source
git clone https://github.com/ethz-adrl/ifopt.git && cd ifopt
mkdir build && cd build
# handles the CMAKE_PREFIX_PATH for the packages installed via vcpkg
cmake --trace-expand -S .. -DCMAKE_TOOLCHAIN_FILE=../../vcpkg/scripts/buildsystems/vcpkg.cmake
make
make install
