![example branch parameter](https://github.com/bsun1220/polytopewalk/actions/workflows/ciwheels.yml/badge.svg?branch=main)

# PolytopeWalk
**PolytopeWalk** is a `C++` library for running MCMC sampling algorithms to generate samples from a uniform distribution over a polytope with a `Python` interface. Existing implementations include the Dikin Walk, John Walk, Vaidya Walk, Ball Walk, Weighted Dikin Walk, and Hit and Run Walk.

## Developer Installation Instructions 
We need to install certain package prerequisites for it to work (listed in each of the operating systems)
- macOS: ``brew install eigen ipopt``
- Windows: ``choco install eigen ipopt``
- Linux: ``yum install -y epel-release eigen3-devel coin-or-Ipopt-devel``

We need to install ifopt using github: 
```
git clone https://github.com/ethz-adrl/ifopt.git && cd ifopt
mkdir build && cd build
cmake ..
make
sudo make install
```

Finally, we can install **PolytopeWalk** as such:
```
git clone https://github.com/bsun1220/polytopewalk
cd polytopewalk
pip install .
```