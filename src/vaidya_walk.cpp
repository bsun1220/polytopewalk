#include "vaidya_walk.hpp"

void VaidyaWalk::generateWeight(VectorXd& x){
    generateDikinHessian(x);
    MatrixXd hess_inv = dhess.inverse();
    MatrixXd slack_inv = slack.cwiseInverse().asDiagonal().toDenseMatrix();

    MatrixXd weights_mat = slack_inv * A * hess_inv * A.transpose() * slack_inv;
    VectorXd wi = weights_mat.diagonal();

    wi.array() += (double)(A.cols()/A.rows());

    weights = wi.asDiagonal().toDenseMatrix();
}

void VaidyaWalk::generateDikinHessian(VectorXd& x){
    generateSlack(x);
    MatrixXd slack_inv = slack.cwiseInverse().asDiagonal().toDenseMatrix();
    dhess = A.transpose() * slack_inv * slack_inv * A;
}

void VaidyaWalk::printType(){
    cout << "Vaidya Walk" << endl;
}