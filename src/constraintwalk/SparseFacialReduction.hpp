#ifndef SPARSE_FR_HPP
#define SPARSE_FR_HPP

#include "SparseLP.hpp"

struct z_res{
    bool found_sol; 
    VectorXd z; 
};

struct fr_res{
    SparseMatrixXd A;
    VectorXd b; 
};

struct res{
    SparseMatrixXd reduced_A;
    VectorXd reduced_b; 
    bool reduced;
};

class SparseFacialReduction {
    public:
        SparseFacialReduction(){}
        res reduce(SparseMatrixXd A, VectorXd b, int k);
        fr_res entireFacialReductionStep(SparseMatrixXd A, VectorXd b, int k);
        z_res findZ(const SparseMatrixXd& A, const VectorXd& b, int k);
    
    protected:
        SparseMatrixXd pickV(const VectorXd& z, int k);
        SparseMatrixXd pickP(const SparseMatrixXd& AV);
};

#endif