#include "SparseDikinLSWalk.hpp"
SparseMatrixXd SparseDikinLSWalk::generateWeight(
    const VectorXd& x, 
    const SparseMatrixXd& A,
    int k
){
    LeverageScore L;

    double d = A.cols() - A.rows();
    double n = k; 
    double q = 2.0 * (1.0 + log(n));
    double alpha = 1.0 - 2.0/q; 

    VectorXd w = VectorXd::Ones(x.rows());
    VectorXd term1 = (alpha) * VectorXd::Ones(x.rows());
    VectorXd errors = ERR * VectorXd::Ones(x.rows());

    for(int i = 0; i < w.rows() - k; i++){
        w(i) = 0;
        term1(i) = 0;
    }
    for(int i = 0; i < MAX_ITER; i++){
        SparseMatrixXd W (x.rows(), x.rows());
        VectorXd term2a = VectorXd::Zero(x.rows());
        for(int j = x.rows() - k; j < x.rows(); j++){
            term2a(j) = (double)alpha/w(j);
            W.coeffRef(j, j) = pow(w(j), alpha * 0.5);
        }
    
        VectorXd term2b = L.generate(A, W, x, ERR, k);
        VectorXd term2 = term2a.cwiseProduct(term2b); 
        VectorXd grad = term1 - term2;
        if (grad.norm() < G_LIM){
            break; 
        }
        w = (w - STEP_SIZE * grad);
        for(int j = x.rows() - k; j < x.rows(); j++){
            w(j) = max(w(j), ERR);
        }
    }
    return SparseMatrixXd(w.asDiagonal());

}

void SparseDikinLSWalk::setDistTerm(int d, int n){
    double q = 2.0 * (1.0 + log(n));
    double term = (1.0 + q) * (1.0 + q * q);
    DIST_TERM = (R * R)/term;
}