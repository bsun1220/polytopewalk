#!/bin/bash

# Clone and bootstrap vcpkg
git clone https://github.com/microsoft/vcpkg
./vcpkg/bootstrap-vcpkg.bat
./vcpkg/vcpkg integrate install

# Install eigen3 and ipopt with vcpkg
./vcpkg/vcpkg install eigen3
# ./vcpkg/vcpkg install coin-or-ipopt

# Test ipopt install from source
svn co https://projects.coin-or.org/svn/Ipopt/stable/3.11 CoinIpopt
cd CoinIpopt/ThirdParty
cd Blas && ./get.Blas && cd ..
cd Lapack && ./get.Lapack && cd ..
cd Metis && ./get.Metis && cd .. (Graph Coloring tool used by e.g. Mumps)
cd Mumps && ./get.Mumps && cd .. (Sparse direct linear solver with permissive license)
cd ..; mkdir build; cd build
../configure --prefix=/usr/local ADD_FFLAGS=-fPIC ADD_CFLAGS=-fPIC ADD_CXXFLAGS=-fPIC
make; sudo make install
pkg-config --libs ipopt

# Install ifopt from source
git clone https://github.com/ethz-adrl/ifopt.git && cd ifopt
mkdir build && cd build
# handles the CMAKE_PREFIX_PATH for the packages installed via vcpkg
cmake --trace-expand -S .. -DCMAKE_TOOLCHAIN_FILE=../../vcpkg/scripts/buildsystems/vcpkg.cmake
make
make install
