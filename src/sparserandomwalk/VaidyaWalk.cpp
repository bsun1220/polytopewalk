#include "VaidyaWalk.hpp"

void VaidyaWalk::generateWeight(const VectorXd& x, const SparseMatrixXd& A, const VectorXd& b){
    generateDikinHessian(x, A, b);

    SimplicialLLT<SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> cholesky;
    cholesky.analyzePattern(dhess);
    cholesky.factorize(dhess);

    SparseMatrixXd I (x.rows(), x.rows());
    for (int i = 0; i < x.rows(); i++){
        I.coeffRef(i, i) = 1;
    }

    SparseMatrixXd hess_inv = cholesky.solve(I);
    SparseMatrixXd slack_inv = slack.cwiseInverse().asDiagonal().toDenseMatrix().sparseView();
    SparseMatrixXd weights_mat = slack_inv * A * hess_inv * A.transpose() * slack_inv;

    VectorXd wi = weights_mat.diagonal();

    wi.array() += (double)(A.cols()/A.rows());

    weights = wi.asDiagonal().toDenseMatrix().sparseView();
}

void VaidyaWalk::generateDikinHessian(const VectorXd& x, const SparseMatrixXd& A, const VectorXd& b){
    generateSlack(x, A, b);
    SparseMatrixXd slack_inv = slack.cwiseInverse().asDiagonal().toDenseMatrix().sparseView();
    dhess = A.transpose() * slack_inv * slack_inv * A;
}

void VaidyaWalk::printType(){
    cout << "Vaidya Walk" << endl;
}