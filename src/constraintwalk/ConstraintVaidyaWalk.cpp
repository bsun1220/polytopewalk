#include "ConstraintVaidyaWalk.hpp"

SparseMatrixXd ConstraintVaidyaWalk::generateG(
    const VectorXd& x, 
    const SparseMatrixXd& A,
    int k
){
    SparseMatrixXd GD = SparseMatrixXd(x.cwiseInverse().cwiseAbs2().asDiagonal());
    SparseMatrixXd GD_inv_sqrt = SparseMatrixXd(x.asDiagonal());
    SparseMatrixXd AGD_inv_sqrt = A * GD_inv_sqrt;
    
    SparseLU<SparseMatrixXd> chol (AGD_inv_sqrt * AGD_inv_sqrt.transpose());
    SparseMatrixXd P = AGD_inv_sqrt.transpose() * chol.solve(AGD_inv_sqrt);
    SparseMatrixXd I = SparseMatrixXd(VectorXd::Ones(A.cols()).asDiagonal());

    SparseMatrixXd M_inv = GD_inv_sqrt * (I - P) * GD_inv_sqrt; 
    SparseMatrixXd slack_inv = SparseMatrixXd(x.cwiseInverse().asDiagonal());

    int n = A.cols();
    int m = A.rows(); 
    SparseMatrixXd res = SparseMatrixXd(slack_inv * M_inv * slack_inv);
    
    VectorXd w = res.diagonal();
    w.array() += double(n - m)/n;

    SparseMatrixXd W = SparseMatrixXd(w.asDiagonal());
    for(int i = 0; i < n - k; i++){
        W.coeffRef(i, i) = ERR; 
    }
    return W;
}

void ConstraintVaidyaWalk::setDistTerm(int d, int n){
    DIST_TERM = (R * R)/sqrt(n * d);
}