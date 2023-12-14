#include "DikinWalk.hpp"

void DikinWalk::setDistTerm(int d, int n){
    
    DIST_TERM = R*R/d;
}

void DikinWalk::generateWeight(const VectorXd& x, const SparseMatrixXd& A, const VectorXd& b){
    int d = b.rows();
    SparseMatrixXd w (d, d);
    for (int i = 0; i < d; i++){
        w.coeffRef(i, i) = 1;
    }
    weights = w;
}

void DikinWalk::printType(){
    cout << "Dikin Walk" << endl;
}