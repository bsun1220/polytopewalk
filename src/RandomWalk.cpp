#include "RandomWalk.hpp"

VectorXd RandomWalk::vectPow(VectorXd& x, float alpha){
    for(int i = 0; i < x.rows(); i++){
        x(i) = pow(x(i), alpha);
    }
    return x;
}

void RandomWalk::initialize(MatrixXd A_p, VectorXd b_p, float r_p){
    A = A_p;
    b = b_p;
    r = r_p;
}

bool RandomWalk::acceptReject(VectorXd& vec, MatrixXd& A, VectorXd&b){
    return ((A * vec) - b).maxCoeff() <= 0;
}

VectorXd RandomWalk::generateGaussianRVNorm(int d){
    VectorXd v(d);
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dis(0.0, 1.0);
    for(int i = 0; i < d; i++){
        v(i) = dis(gen);
    }
    return v/v.norm();
}

MatrixXd RandomWalk::generateCompleteWalk(int num_steps, VectorXd& x){
    cout << "oops" << endl;
    return MatrixXd::Zero(1,1);
}