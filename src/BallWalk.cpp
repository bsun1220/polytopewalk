#include "BallWalk.hpp"


MatrixXd BallWalk::generateCompleteWalk(const int num_steps, VectorXd& x){
    int n = x.rows(); 
    MatrixXd results = MatrixXd::Zero(num_steps, n);
    for (int i = 0; i < num_steps; i++){
        VectorXd new_x = generateGaussianRVNorm(n) * r + x;
        if (acceptReject(new_x, A, b)){
            x = new_x;
        }
        results.row(i) = x; 
    }
    return results;
}

void BallWalk::printType(){
    cout << "Ball Walk" << endl;
}