#!/bin/bash

# Clone and bootstrap vcpkg
git clone https://github.com/microsoft/vcpkg
./vcpkg/bootstrap-vcpkg.bat
./vcpkg/vcpkg integrate install

# Install necessary libraries with vcpkg
./vcpkg/vcpkg install eigen3
./vcpkg/vcpkg install ipopt
