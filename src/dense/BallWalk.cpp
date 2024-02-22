#include "BallWalk.hpp"


MatrixXd BallWalk::generateCompleteWalk(const int num_steps, VectorXd& x, const MatrixXd& A, const VectorXd& b){
    int n = x.rows(); 
    int d = A.cols();
    MatrixXd results = MatrixXd::Zero(num_steps, n);
    for (int i = 0; i < num_steps; i++){
        VectorXd new_x = generateGaussianRVNorm(n) * R/sqrt(d) + x;
        if (inPolytope(new_x, A, b)){
            x = new_x;
        }
        results.row(i) = x; 
    }
    return results;
}

void BallWalk::printType(){
    cout << "Ball Walk" << endl;
}