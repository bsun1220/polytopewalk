#!/bin/bash

# Clone and bootstrap vcpkg
git clone https://github.com/microsoft/vcpkg
./vcpkg/bootstrap-vcpkg.bat
./vcpkg/vcpkg integrate install

# Install eigen3 and ipopt with vcpkg
./vcpkg/vcpkg install eigen3
./vcpkg/vcpkg install coin-or-ipopt

# get FindIPOPT_DIR from casadi, which is better written
git clone --depth 1 --branch 3.6.5 https://github.com/casadi/casadi.git

# install ifopt from source
git clone https://github.com/ethz-adrl/ifopt.git
cd ifopt
# move FindIPOPT.cmake around
mv ifopt_ipopt/cmake/FindIPOPT.cmake ifopt_ipopt/cmake/FindIPOPT.cmakeold
cp ../casadi/cmake/FindIPOPT.cmake ifopt_ipopt/cmake/
cp ../casadi/cmake/canonicalize_paths.cmake ifopt_ipopt/cmake/
mkdir build
cd build
cmake .. \
  -DCMAKE_TOOLCHAIN_FILE="../../vcpkg/scripts/buildsystems/vcpkg.cmake" \
  -DIPOPT_LIBRARIES="../../vcpkg/installed/x64-windows/lib/ipopt.lib" \
  -DIPOPT_INCLUDE_DIRS="../../vcpkg/installed/x64-windows/include/coin-or" \
  -G "Unix Makefiles"

make
make install
cd ..
cd ..

export CMAKE_TOOLCHAIN_FILE="../../vcpkg/scripts/buildsystems/vcpkg.cmake"

