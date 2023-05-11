#include "DikinWalk.hpp"

void DikinWalk::generateWeight(VectorXd& x){
    int d = b.rows();
    weights = MatrixXd::Identity(d, d);
}

void DikinWalk::printType(){
    cout << "Dikin Walk" << endl;
}