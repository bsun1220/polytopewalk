#include "VaidyaWalk.hpp"

void VaidyaWalk::generateWeight(const VectorXd& x, const SparseMatrixXd& A, const VectorXd& b){
    generateDikinHessian(x, A, b);

    SparseMatrixXd I (x.rows(), x.rows());
    for (int i = 0; i < x.rows(); i++){
        I.coeffRef(i, i) = 1;
    }

    SparseLU<SparseMatrixXd> chol (dhess);

    SparseMatrixXd hess_inv = chol.solve(I);
    SparseMatrixXd slack_inv = SparseMatrixXd(slack.cwiseInverse().asDiagonal());

    SparseMatrixXd weights_mat = slack_inv * A * hess_inv * A.transpose() * slack_inv;

    VectorXd wi = weights_mat.diagonal();

    wi.array() += (double)(A.cols()/A.rows());

    weights = SparseMatrixXd(wi.asDiagonal());

}

void VaidyaWalk::generateDikinHessian(const VectorXd& x, const SparseMatrixXd& A, const VectorXd& b){
    generateSlack(x, A, b);
    SparseMatrixXd slack_inv = SparseMatrixXd(slack.cwiseInverse().asDiagonal());
    dhess = A.transpose() * slack_inv * slack_inv * A;
}

void VaidyaWalk::printType(){
    cout << "Vaidya Walk" << endl;
}