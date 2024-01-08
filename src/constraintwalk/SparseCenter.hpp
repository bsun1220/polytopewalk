#ifndef SPARSE_CENTER_HPP
#define SPARSE_CENTER_HPP

#include "Common.hpp"
#include "SparseLP.hpp"

class SparseCenter {
    public:
        SparseCenter(){}
        VectorXd getInitialPoint(SparseMatrixXd A, VectorXd b, int k);
};

#endif