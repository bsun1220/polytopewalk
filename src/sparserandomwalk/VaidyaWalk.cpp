#include "VaidyaWalk.hpp"

void VaidyaWalk::generateWeight(const VectorXd& x, const SparseMatrixXd& A, const VectorXd& b){
    generateDikinHessian(x, A, b);


    SparseMatrixXd A_transpose = A.transpose();

    SimplicialCholesky<SparseMatrixXd> chol (dhess);
    SparseMatrixXd hess_inv = chol.solve(A_transpose);
    SparseMatrixXd slack_inv = SparseMatrixXd(slack.cwiseInverse().asDiagonal());

    SparseMatrixXd weights_mat = slack_inv * A * hess_inv * slack_inv;

    VectorXd wi = weights_mat.diagonal();

    wi.array() += (double)(A.cols()/A.rows());

    weights = SparseMatrixXd(wi.asDiagonal());

}

void VaidyaWalk::generateDikinHessian(const VectorXd& x, const SparseMatrixXd& A, const VectorXd& b){
    generateSlack(x, A, b);
    VectorXd slack_inv_v = slack.cwiseInverse();
    SparseMatrixXd slack_inv = SparseMatrixXd(slack.cwiseInverse().asDiagonal());
    dhess = A.transpose() * slack_inv * slack_inv * A;
}

void VaidyaWalk::printType(){
    cout << "Vaidya Walk" << endl;
}