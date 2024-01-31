#include "ConstraintVaidyaWalk.hpp"
#include "LeverageScore.hpp"

SparseMatrixXd ConstraintVaidyaWalk::generateG(
    const VectorXd& x, 
    const SparseMatrixXd& A,
    int k
){
    LeverageScore L;
    SparseMatrixXd W (x.rows(), x.rows());
    for(int i = 0; i < x.rows(); i++){
        W.coeffRef(i, i) = 1;
    }
    VectorXd weights = L.generate(A, W, x);
    for(int i = 0; i < weights.rows() - k; i++){
        weights(i) = 0; 
    }
    for (int i = weights.rows() - k; i < weights.rows(); i++){
        weights(i) += ((double)(A.cols() - A.rows())/k);
    }
    return SparseMatrixXd(weights.asDiagonal());
}

void ConstraintVaidyaWalk::setDistTerm(int d, int n){
    DIST_TERM = (R * R)/sqrt(n * d);
}