#!/bin/bash

# install prerequisite
pacman -Syu
# pacman -S binutils diffutils git grep make patch pkg-config
pacman -S git
pacman -S mingw-w64-x86_64-lapack
pacman -S mingw-w64-x86_64-metis
pacman -S mingw-w64-x86_64-eigen3

echo "OK pause"


