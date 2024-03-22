#include "SparseJohnWalk.hpp"
#include "LeverageScore.hpp"

SparseMatrixXd SparseJohnWalk::generateWeight(
    const VectorXd& x, 
    const SparseMatrixXd& A,
    int k
){
    double d = A.cols() - A.rows();
    double n = k; 
    double alpha = 1 - 1/(log2(2.0 * n / d));
    double beta = (double)d / (2.0 * n);

    LeverageScore L;
    VectorXd w = VectorXd::Ones(x.rows());
    VectorXd beta_ones = beta * VectorXd::Ones(x.rows());
    for(int i = 0; i < w.rows() - k; i++){
        w(i) = 0;
        beta_ones.coeffRef(i) = 0;
    }
    VectorXd next_weight = w; 

    for(int i = 0; i < MAX_ITER; i++){
        w = next_weight;
        SparseMatrixXd W (w.rows(), w.rows());
        for(int j = x.rows() - k; j < x.rows(); j++){
            W.coeffRef(j, j) = pow(w(j), alpha * 0.5);
        }
        VectorXd score  = L.generate(A, W, x, ERR, k);
        next_weight = 0.5 * (w + score + beta_ones).cwiseMax(beta_ones);

        if ((w - next_weight).cwiseAbs().maxCoeff() < LIM){
            break;
        }
    }

    return SparseMatrixXd(w.asDiagonal());

}

void SparseJohnWalk::setDistTerm(int d, int n){
    DIST_TERM = (R * R)/pow(d, 1.5);
}