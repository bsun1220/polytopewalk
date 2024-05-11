#!/bin/bash

# install prerequisite
# pacman -Syu
# pacman -S binutils diffutils git grep make patch pkg-config
# pacman -S mingw-w64-x86_64-lapack
# pacman -S mingw-w64-x86_64-metis
# pacman -S mingw-w64-x86_64-eigen3

# pacman -Qs eigen3
# pacman -Ql eigen3 | grep bin/

ls "C:\msys64"
ls "C:\msys64\mingw64"
ls "C:\msys64\mingw64\bin"


echo "OK pause"
export PATH="/c/Program Files/Git/bin:$PATH"

git clone https://github.com/ethz-adrl/ifopt.git
cd ifopt
mkdir build
cd build
cmake -S ..
make install



