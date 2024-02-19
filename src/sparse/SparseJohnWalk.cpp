#include "SparseJohnWalk.hpp"
#include "LeverageScore.hpp"

SparseMatrixXd SparseJohnWalk::generateWeight(
    const VectorXd& x, 
    const SparseMatrixXd& A,
    int k
){
    LeverageScore L;
    VectorXd term1 = VectorXd::Ones(x.rows());
    VectorXd w = VectorXd::Ones(x.rows());
    for(int i = 0; i < w.rows() - k; i++){
        w(i) = 0;
        term1(i) = 0;
    }

    double d = A.cols() - A.rows();
    double n = k; 
    double alpha = 1 - 1/(log2(2.0 * n / d));
    double beta = (double)d / (2.0 * n);

    for(int i = 0; i < MAX_ITER; i++){
        SparseMatrixXd W (w.rows(), w.rows());
        VectorXd term2a = VectorXd::Zero(w.rows());
        for(int j = x.rows() - k; j < x.rows(); j++){
            term2a(j) = 1.0/w(j);
            W.coeffRef(j, j) = pow(w(j), alpha * 0.5);
        }

        VectorXd term2b = L.generate(A, W, x, ERR, k);
        VectorXd term2 = term2a.cwiseProduct(term2b); 
        VectorXd term3 = VectorXd::Zero(w.rows());
        for(int j = w.rows() - k; j < w.rows(); j++){
            term3(j) = beta * 1.0/w(j);
        }

        VectorXd grad = term1 - term2 - term3;
        if (grad.norm() < G_LIM){
            break; 
        }
        w = w - STEP_SIZE * grad; 
        for(int j = x.rows() - k; j < x.rows(); j++){
            w(j) = max(w(j), beta);
        }

    }
    return SparseMatrixXd(w.asDiagonal());

}

void SparseJohnWalk::setDistTerm(int d, int n){
    DIST_TERM = (R * R)/pow(d, 1.5);
}