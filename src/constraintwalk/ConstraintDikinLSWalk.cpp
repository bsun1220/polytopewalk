#include "ConstraintDikinLSWalk.hpp"
SparseMatrixXd ConstraintDikinLSWalk::generateG(
    const VectorXd& x, 
    const SparseMatrixXd& A,
    int k
){
    LeverageScore L;
    VectorXd w = VectorXd::Ones(x.rows());

    for(int i = 0; i < w.rows() - k; i++){
        w(i) = 0; 
    }

    double d = A.cols() - A.rows();
    double n = k; 
    double q = 2.0 * (1.0 + log(n));
    double alpha = 1.0 - 2.0/q; 
    VectorXd term1 = (0.5 - 1.0/q) * VectorXd::Ones(x.rows());
    double beta = (double)d / (2.0 * n);

    for(int i = 0; i < MAX_ITER; i++){
        SparseMatrixXd W (w.rows(), w.rows());
        for(int j = 0; j < w.rows(); j++){
            W.coeffRef(j, j) = pow(w(j), alpha);
        }
        VectorXd term2a = alpha * w.array().pow(alpha - 1);
        VectorXd term2b = L.generate(A, W, x);
        VectorXd term2 = term2a.asDiagonal() * term2b; 
        VectorXd term3 = beta * w.cwiseInverse();
        VectorXd grad = term1 - term2 - term3;

        for(int i = 0; i < w.rows() - k; i++){
            grad(i) = 0; 
        }

        w = w - STEP_SIZE * grad; 

        if (grad.norm() < G_LIM){
            break; 
        }

    }
    return SparseMatrixXd(w.asDiagonal());

}

void ConstraintDikinLSWalk::setDistTerm(int d, int n){
    double q = 2.0 * (1.0 + log(n));
    double term = (1.0 + q) * (1.0 + q * q);
    DIST_TERM = (R * R)/term;
}