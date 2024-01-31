#include "ConstraintJohnWalk.hpp"
#include "LeverageScore.hpp"

SparseMatrixXd ConstraintJohnWalk::generateG(
    const VectorXd& x, 
    const SparseMatrixXd& A,
    int k
){
    LeverageScore L;
    VectorXd term1 = VectorXd::Ones(x.rows());
    VectorXd w = VectorXd::Ones(x.rows());
    for(int i = 0; i < w.rows() - k; i++){
        w(i) = 0; 
    }
    double d = A.cols() - A.rows();
    double n = k; 
    double alpha = 1 - 1/(log2(2.0 * n / d));
    double beta = (double)d / (2.0 * n);

    for(int i = 0; i < MAX_ITER; i++){
        SparseMatrixXd W (w.rows(), w.rows());
        for(int j = 0; j < w.rows(); j++){
            W.coeffRef(j, j) = pow(w(j), alpha);
        }
        VectorXd term2a = w.array().pow(alpha - 1);
        VectorXd term2b = L.generate(A, W, x);
        VectorXd term2 = term2a.asDiagonal() * term2b; 
        VectorXd term3 = beta * w.cwiseInverse();
        VectorXd grad = term1 - term2 - term3;

        for(int i = 0; i < w.rows() - k; i++){
            grad(i) = 0; 
        }
        if (grad.norm() < G_LIM){
            break; 
        }
        w = w - STEP_SIZE * grad; 

    }
    return SparseMatrixXd(w.asDiagonal());

}

void ConstraintJohnWalk::setDistTerm(int d, int n){
    DIST_TERM = (R * R)/pow(d, 1.5);
}