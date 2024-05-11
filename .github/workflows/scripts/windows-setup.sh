#!/bin/bash

# Clone and bootstrap vcpkg
git clone https://github.com/microsoft/vcpkg
./vcpkg/bootstrap-vcpkg.bat
./vcpkg/vcpkg integrate install

cd vcpkg
# Set the VCPKG_ROOT environment variable to the current directory
export VCPKG_ROOT=$(pwd)
# Display VCPKG_ROOT to verify it is set correctly
echo "VCPKG_ROOT is set to: $VCPKG_ROOT"
cd ..


# Install eigen3 and ipopt with vcpkg
./vcpkg/vcpkg install eigen3
./vcpkg/vcpkg install coin-or-ipopt

# export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:./vcpkg/installed/x64-windows/lib/ipopt.lib"

echo `ls ./vcpkg/`
echo `ls ./vcpkg/installed`
echo `ls ./vcpkg/installed/x64-windows`
echo `ls ./vcpkg/installed/x64-windows/lib`

# get FindIPOPT_DIR.cmake from casadi, which is better written
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
  -DCMAKE_TOOLCHAIN_FILE="$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" \
  -DIPOPT_LIBRARIES="$VCPKG_ROOT/installed/x64-windows/lib/ipopt.lib" \
  -DIPOPT_INCLUDE_DIRS="$VCPKG_ROOT/installed/x64-windows/include/coin-or/" \
  -G "Unix Makefiles"

make
make install
cd ..
cd ..

