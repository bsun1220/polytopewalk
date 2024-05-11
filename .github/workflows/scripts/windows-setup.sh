#!/bin/bash

# install prerequisite
pacman -Syu
# pacman -S binutils diffutils git grep make patch pkg-config
pacman -S mingw-w64-x86_64-lapack
pacman -S mingw-w64-x86_64-metis
pacman -S mingw-w64-x86_64-eigen3

pacman -Qs eigen3
pacman -Ql eigen3 | grep bin/


echo "OK pause"

C:\Program Files\Git\bin\git.exe clone https://github.com/ethz-adrl/ifopt.git
cd ifopt
mkdir build
cd build
cmake -S ..
make install



