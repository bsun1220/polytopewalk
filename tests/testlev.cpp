#include <sparserandomwalk/PolytopeWalk.hpp>
#include <set>

int main(){


    MatrixXd At (4,2);
    At << 1, 0, -1, 0, 0, 1, 0, -1;
    SparseMatrixXd A = At.sparseView();
    VectorXd b (4);
    b << 1,1,1,1;
    VectorXd x (2);
    x << 0.4, 0.3;

    VectorXd slack = (b - A * x);

    SparseMatrixXd slack_inv = SparseMatrixXd(slack.cwiseInverse().asDiagonal());
    SparseMatrixXd hess = A.transpose() * slack_inv * slack_inv * A;
    SparseMatrixXd I (x.rows(), x.rows());
    for (int i = 0; i < x.rows(); i++){
        I.coeffRef(i, i) = 1;
    }
    
    SparseLU<SparseMatrixXd> chol (hess);
    SparseMatrixXd hess_inv = chol.solve(I);
    SparseMatrixXd weights = (slack_inv * A * hess_inv * A.transpose() * slack_inv);
    VectorXd lev = weights.diagonal();


    SparseMatrixXd Ia (A.rows(), A.rows());
    for (int i = 0; i < A.rows(); i++){
        Ia.coeffRef(i, i) = 1;
    }

    LeverageScore L;

    VectorXd pred = L.generate(A, Ia, b, x);

    cout << "True" << endl;
    cout << lev << endl;
    cout << "Predicted" << endl;
    cout << pred << endl;

}