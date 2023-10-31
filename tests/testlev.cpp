#include <sparserandomwalk/PolytopeWalk.hpp>
#include <set>

int main(){

    /*
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
    VectorXd lev = weights.diagonal();*/

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    int dim = 4;
    MatrixXd dhess (dim, dim);
    for(int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            dhess.coeffRef(i, j) = dis(gen);
        }
    }
    MatrixXd dhess_transpose = dhess.transpose();
    dhess = (dhess + dhess_transpose)/2;
    dhess = dhess + dim * MatrixXd::Identity(dim, dim);

    MatrixXd hess_inv = dhess.inverse();

    SparseMatrixXd hess = dhess.sparseView();
    SimplicialLDLT<SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> cholesky;
    cholesky.analyzePattern(hess);
    cholesky.factorize(hess);

    SparseMatrixXd L0 = cholesky.matrixL();
    SparseMatrixXd U0 = cholesky.matrixU();
    VectorXd D = cholesky.vectorD();
    SparseMatrixXd D_sqrt (L0.rows(), U0.cols());


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
    map<int, vector<pair<int, int>>>col_major;


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

    cout << "True" << endl;
    cout << MatrixXd(hess_inv) << endl;
    cout << "Predicted" << endl;
    cout << MatrixXd(inv) << endl;

}