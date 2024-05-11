#!/bin/bash

# install prerequisite
# pacman -Syu
# pacman -S binutils diffutils git grep make patch pkg-config
# pacman -S mingw-w64-x86_64-lapack
# pacman -S mingw-w64-x86_64-metis
# pacman -S mingw-w64-x86_64-eigen3

# pacman -Qs eigen3
# pacman -Ql eigen3 | grep bin/

echo "Let's start windows-setup"
export PATH="/c/msys64/mingw64/bin:/c/Program Files/Git/bin:$PATH"
# echo "PATH=$PATH:/c/msys64/mingw64/bin" >> $GITHUB_ENV

# install ipopt via https://coin-or.github.io/Ipopt/INSTALL.html
# install Mumps
git clone https://github.com/coin-or-tools/ThirdParty-Mumps.git
cd ThirdParty-Mumps
./get.Mumps
./configure
make
make install
cd ..

echo `pwd`

# install ipopt from source
git clone https://github.com/coin-or/Ipopt.git
cd Ipopt
mkdir build
cd build
../configure
make
# make test
make install
export IPOPT_DIR=`pwd`
cd ..

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/mingw64/lib:/c/msys64/mingw64/lib"

echo $IPOPT_DIR
# echo "list d"
# ls "/d/a/_temp/msys64/mingw64/include"
echo "list c222222"
ls "/c/msys64/mingw64/lib"


# install ifopt from source
git clone https://github.com/ethz-adrl/ifopt.git
cd ifopt
mkdir build
cd build
export CMAKE_PREFIX_PATH="/mingw64:/c/msys64/mingw64:$CMAKE_PREFIX_PATH"
echo $CMAKE_PREFIX_PATH
cmake ..
make install
cd ..



