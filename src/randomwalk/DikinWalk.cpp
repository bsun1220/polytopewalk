#include "DikinWalk.hpp"

void DikinWalk::generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    int d = b.rows();
    weights = MatrixXd::Identity(d, d);
}

void DikinWalk::printType(){
    cout << "Dikin Walk" << endl;
}