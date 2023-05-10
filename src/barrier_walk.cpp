#include "barrier_walk.hpp"

void BarrierWalk::initialize(MatrixXd A_p, VectorXd b_p, float r){
    A = A_p;
    b = b_p;
    r = r; 
}

void BarrierWalk::setTs(float a){
    term_sample = a; 
}

void BarrierWalk::setTd(float b){
    term_density = b; 
}

VectorXd BarrierWalk::generateGaussianRV(int d){
    VectorXd v(d);
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dis(0.0, 1.0);
    for(int i = 0; i < d; i++){
        v(i) = dis(gen);
    }
    return v;
}

void BarrierWalk::generateSlack(VectorXd& x){
    slack = (b - (A * x));
}

float BarrierWalk::localNorm(VectorXd v, MatrixXd& m){
    return ((v.transpose() * m) * v)(0);
}

void BarrierWalk::generateWeight(VectorXd& x){
    int d = b.rows();
    weights = VectorXd::Zero(d).asDiagonal().toDenseMatrix();
}

void BarrierWalk::generateHessian(VectorXd& x){
    generateWeight(x);
    generateSlack(x);
    MatrixXd slack_inv = slack.cwiseInverse().asDiagonal().toDenseMatrix();
    hess = A.transpose() * slack_inv * weights * slack_inv * A;
}

float BarrierWalk::generateProposalDensity(VectorXd& x, VectorXd& z){
    generateHessian(x);
    VectorXd d = generateGaussianRV(x.rows());
    return sqrt(hess.determinant()) * exp(term_density * localNorm(x - z, hess));
}

void BarrierWalk::generateSample(VectorXd& x){
    generateHessian(x);
    MatrixXd matrix = hess.inverse().sqrt();
    VectorXd direction = generateGaussianRV(x.rows());
    z = x + term_sample * (matrix * direction);
}

void BarrierWalk::printType(){
    cout << "Generic Barrier" << endl;
}

MatrixXd BarrierWalk::generateCompleteWalk(int num_steps, VectorXd& x){
    MatrixXd results = MatrixXd::Zero(num_steps, A.cols());
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);
    float one = 1.0;
    for(int i = 0; i < num_steps; i++){
        generateSample(x);
        if(acceptReject(z, A, b)){
            float g_x_z = generateProposalDensity(x, z);
            float g_z_x = generateProposalDensity(z, x);
            float alpha = min(one, g_z_x/g_x_z);
            float val = dis(gen);
            x = val < alpha ? z : x;
        }
        results.row(i) = x.transpose();
    }
    return results;
}