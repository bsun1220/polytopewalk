#include <sparserandomwalk/PolytopeWalk.hpp>

int main(){
    MatrixXd At (4,2);
    At << 1, 0, -1, 0, 0, 1, 0, -1;
    SparseMatrixXd A = At.sparseView();
    VectorXd b (4);
    b << 1,1,1,1;
    VectorXd x (2);
    x << 0.5, 0.5;

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

    SimplicialLLT<SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> cholesky;
    cholesky.analyzePattern(hess);
    cholesky.factorize(hess);
    SparseMatrixXd L = cholesky.matrixL();

    vector<pair<int, int>> lst; 

    for(int k = 0; k < L.outerSize(); k++){
        for(SparseMatrixXd::InnerIterator it(L, k); it; ++it){
            lst.push_back(make_pair(it.row(), it.col()));
        }
    }

    SparseMatrixXd S (hess.rows(), hess.cols());
    for(pair<int, int> s : lst){
        int row = s.first;
        int col = s.second; 
        S.insert(row, col) = hess.coeff(row, col);

    }
    cout << MatrixXd(hess) << endl;
    cout << MatrixXd(S) << endl;

    cout << "****"<<endl;
    cout << lev << endl;

}