#!/bin/bash

# install prerequisite
pacman -Syu
# pacman -S binutils diffutils git grep make patch pkg-config
pacman -S git
pacman -S gcc fortran
pacman -S lapack metis


