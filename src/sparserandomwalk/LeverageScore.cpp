#include "LeverageScore.hpp"

VectorXd LeverageScore::generate(const SparseMatrixXd& A, const SparseMatrixXd& W, const VectorXd& b, const VectorXd& x){

    VectorXd slack = (b - (A * x));

    SparseMatrixXd slack_inv = SparseMatrixXd(slack.cwiseInverse().asDiagonal());

    SparseMatrixXd hess = A.transpose() * slack_inv * W * slack_inv * A;
    SimplicialLDLT<SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> cholesky;
    cholesky.analyzePattern(hess);
    cholesky.factorize(hess);
    SparseMatrixXd L0 = cholesky.matrixL();
    VectorXd D = cholesky.vectorD();
    SparseMatrixXd D_sqrt (L0.rows(), L0.cols());

    for(int i = 0; i < L0.rows(); i++){
        D_sqrt.coeffRef(i, i) = sqrt(D(i));
    }
    SparseMatrixXd L = L0 * D_sqrt;
    map<int, vector<int>> sparsity_ref; 
    // row : list of col indices
    for(int k = 0; k < L.outerSize(); k++){
        for(SparseMatrixXd::InnerIterator it(L, k); it; ++it){
            sparsity_ref[it.col()].push_back(it.row());
        }
    }
    SparseMatrixXd inv(L.rows(), L.rows());
    map<int, vector<pair<int, int>>> row_major;
    map<int, vector<pair<int, int>>> col_major;

    for(int i = L.rows() - 1; i >= 0; i--){
        for(int v = sparsity_ref[i].size() - 1; v >= 0; v--){
            int j = sparsity_ref[i][v];
            double z = (i == j) ? (double)1/D(i) : 0;

            for(auto pair : col_major[j]){
                int first = pair.first;
                int second = pair.second;
                double term = inv.coeff(first, second);
                int k = min(first, second);

                z -= term * L0.coeff(k, i);
            }

            for(auto pair : row_major[j]){
                int first = pair.first;
                int second = pair.second;
                double term = inv.coeff(first, second);
                int k = max(first, second);
                z -= term * L0.coeff(k, i);
            }

            col_major[j].push_back(make_pair(i, j));
            if (i != j){
                row_major[i].push_back(make_pair(i, j));
            }

            inv.coeffRef(i, j) = z;
        }
    }

    SparseMatrixXd res = slack_inv * A * inv * A.transpose() * slack_inv;
    return res.diagonal();
}
