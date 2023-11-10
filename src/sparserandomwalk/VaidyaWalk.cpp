#include "VaidyaWalk.hpp"

void VaidyaWalk::generateWeight(const VectorXd& x, const SparseMatrixXd& A, const VectorXd& b){

    SparseMatrixXd I(A.rows(), A.rows());
    for(int i = 0; i < A.rows(); i++){
        I.coeffRef(i, i) = 1;
    }
    LeverageScore L;
    VectorXd wi = L.generate(A, I, b, x);
    wi.array() += (double)(A.cols()/A.rows());
    weights = SparseMatrixXd(wi.asDiagonal());
}

void VaidyaWalk::printType(){
    cout << "Vaidya Walk" << endl;
}