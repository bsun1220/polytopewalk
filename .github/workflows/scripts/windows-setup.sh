#!/bin/bash

echo "Let's start windows-setup"
# export PATH="/c/msys64/mingw64/bin:/c/Program Files/Git/bin:$PATH"
export PATH="/mingw64/bin:$PATH"

# install ipopt via https://coin-or.github.io/Ipopt/INSTALL.html
# install Mumps
git clone https://github.com/coin-or-tools/ThirdParty-Mumps.git
cd ThirdParty-Mumps
./get.Mumps
./configure --prefix=/mingw64/
make
make install
cd ..

# install ipopt from source
git clone https://github.com/coin-or/Ipopt.git
cd Ipopt
mkdir build
cd build
../configure --prefix=/mingw64/
make
# make test
make install
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
cmake .. -DCMAKE_VERBOSE_MAKEFILE=ON \
  -DCMAKE_INSTALL_PREFIX="/mingw64/" \
  -DCMAKE_PREFIX_PATH="/mingw64/lib" \
  -DIPOPT_LIBRARIES="/mingw64/lib/libipopt.dll.a" \
  -DIPOPT_INCLUDE_DIRS="/mingw64/include/coin-or" \
  -G "Unix Makefiles"

make VERBOSE=1
make install
cd ..
cd ..

# ls /c/msys64/mingw64/share
# ls /c/msys64/mingw64/share/eigen3
# ls /c/msys64/mingw64/share/ifopt
# ls /mingw64/share
# ls /c/msys64/mingw64/bin
# ls /mingw64/bin

cmake_predix_path=$(cygpath -w /mingw64/share)
echo "CMAKE_PREFIX_PATH=$cmake_predix_path" >> $GITHUB_ENV
echo $cmake_predix_path
eigen_dir=$(cygpath -w /mingw64/share/eigen3/cmake)
echo "Eigen3_DIR=$eigen_dir" >> $GITHUB_ENV
ifopt_dir=$(cygpath -w /mingw64/share/ifopt/cmake)
echo $ifopt_dir
echo "ifopt_DIR=$ifopt_dir" >> $GITHUB_ENV
# # get CXX compiler location
# cc=$(cygpath -w /mingw64/bin/cc.exe)
# echo "CC=$cc" >> $GITHUB_ENV
# cxx=$(cygpath -w /mingw64/bin/c++.exe)
# echo "CXX=$cxx" >> $GITHUB_ENV
