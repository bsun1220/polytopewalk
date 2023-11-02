#include "VaidyaWalk.hpp"

void VaidyaWalk::generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    generateDikinHessian(x, A, b);
    FullPivLU<MatrixXd> lu(dhess);
    MatrixXd I = Eigen::MatrixXd::Identity(x.rows(),x.rows());
    MatrixXd hess_inv = lu.solve(I);

    MatrixXd slack_inv = slack.cwiseInverse().asDiagonal().toDenseMatrix();

    MatrixXd weights_mat = slack_inv * A * hess_inv * A.transpose() * slack_inv;
    VectorXd wi = weights_mat.diagonal();

    wi.array() += (double)(A.cols()/A.rows());

    weights = wi.asDiagonal().toDenseMatrix();
}

void VaidyaWalk::generateDikinHessian(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    generateSlack(x, A, b);
    MatrixXd slack_inv = slack.cwiseInverse().asDiagonal().toDenseMatrix();
    dhess = A.transpose() * slack_inv * slack_inv * A;
}

void VaidyaWalk::printType(){
    cout << "Vaidya Walk" << endl;
}