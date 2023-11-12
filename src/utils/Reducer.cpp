#include "Reducer.hpp"

MatrixXd Reducer::makeFullRank(const MatrixXd& mat){
    if(mat.rows() == mat.cols()){
        return mat; 
    }

    
    HouseholderQR <MatrixXd> qr(mat.cols(), mat.rows());
    qr.compute(mat.transpose());
    MatrixXd q = qr.householderQ();
    MatrixXd r = qr.matrixQR().triangularView<Eigen::Upper>();

    MatrixXd iden = MatrixXd::Identity(q.rows(), q.rows());
    MatrixXd new_r (q.rows(), q.rows());
    for(int i = 0; i < r.cols(); i++){
        new_r.col(i) = r.col(i);
    }
    for(int i = r.cols(); i < q.rows(); i++){
        new_r.col(i) = iden.col(i);
    }

    MatrixXd ans = (q * new_r).transpose();
    
    const double multiplier = std::pow(10.0, 10);
    for(int i = 0; i < ans.rows(); i++){
        for(int j = 0; j < ans.cols(); j++){
            ans.coeffRef(i,j) = round(ans.coeffRef(i,j) * multiplier) / multiplier;
        }
    }
    return ans;
}

problem_result Reducer::reduce( MatrixXd A, VectorXd b){
    problem_result pf;
    cout << "oops" << endl;
    return pf;
}