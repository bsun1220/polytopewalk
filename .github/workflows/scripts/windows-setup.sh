#!/bin/bash

# install prerequisite
pacman -S binutils diffutils git grep make patch pkg-config
pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran
pacman -S mingw-w64-x86_64-lapack mingw-w64-x86_64-metis


