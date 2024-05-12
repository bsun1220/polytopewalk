#!/bin/bash

echo "Let's start windows-setup"
export PATH="/c/msys64/mingw64/bin:/c/Program Files/Git/bin:$PATH"
echo "PATH=$PATH:/c/msys64/mingw64/bin" >> $GITHUB_ENV

# install ipopt via https://coin-or.github.io/Ipopt/INSTALL.html
# install Mumps
git clone https://github.com/coin-or-tools/ThirdParty-Mumps.git
cd ThirdParty-Mumps
./get.Mumps
./configure --prefix=/c/msys64/mingw64/
make
make install
cd ..

# install ipopt from source
git clone https://github.com/coin-or/Ipopt.git
cd Ipopt
mkdir build
cd build
../configure --prefix=/c/msys64/mingw64/
make
# make test
make install
# export IPOPT_DIR=`pwd`
cd ..
cd ..

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
cmake .. -DIPOPT_LIBRARIES="/c/msys64/mingw64/lib/libipopt.dll.a" -DIPOPT_INCLUDE_DIRS="/c/msys64/mingw64/include/coin-or" -G "MinGW Makefiles"

make
make install
cd ..
cd ..


