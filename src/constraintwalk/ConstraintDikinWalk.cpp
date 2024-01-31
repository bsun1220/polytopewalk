#include "ConstraintDikinWalk.hpp"

SparseMatrixXd ConstraintDikinWalk::generateG(
    const VectorXd& x, 
    const SparseMatrixXd& A,
    int k
){
    SparseMatrixXd G (A.cols(), A.cols());
    int n = A.cols(); 

    for(int i = 0; i < n - k; i++){
        G.coeffRef(i, i) = 0;
    }

    for(int i = n - k; i < n; i++){
        G.coeffRef(i, i) = 1/(x[i] * x[i]);
    }
    return G; 
}

void ConstraintDikinWalk::setDistTerm(int d, int n){
    DIST_TERM = (R * R)/d; 
}
