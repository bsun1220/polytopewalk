#include "ConstraintWalk.hpp"

VectorXd ConstraintWalk::generateGaussianRV(int d){
    VectorXd v(d);
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dis(0.0, 1.0);
    for(int i = 0; i < d; i++){
        v(i) = dis(gen);
    }
    return v;
}

MatrixXd ConstraintWalk::generateCompleteWalk(
    const int num_steps,
    const VectorXd& init, 
    const SparseMatrixXd& A,
    const VectorXd& b, 
    int k
){

    cout << "Oops" << endl;
    return MatrixXd::Zero(1,1);

}

bool ConstraintWalk::inPolytope(
    const VectorXd&z, 
    int k
){
    return z.tail(k).minCoeff() >= 0; 
}

bool ConstraintWalk::inApproxPolytope(
    const VectorXd&z, 
    int k,
    double err
){
    return z.tail(k).minCoeff() >= -err;
}