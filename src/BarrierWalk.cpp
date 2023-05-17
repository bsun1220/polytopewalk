#include "BarrierWalk.hpp"


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

void BarrierWalk::generateSlack(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    slack = (b - (A * x));
}

float BarrierWalk::localNorm(VectorXd v, const MatrixXd& m){
    return ((v.transpose() * m) * v)(0);
}

void BarrierWalk::generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    int d = b.rows();
    weights = VectorXd::Zero(d).asDiagonal().toDenseMatrix();
}

void BarrierWalk::generateHessian(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    generateWeight(x, A, b);
    generateSlack(x, A, b);
    MatrixXd slack_inv = slack.cwiseInverse().asDiagonal().toDenseMatrix();
    hess = A.transpose() * slack_inv * weights * slack_inv * A;
}

float BarrierWalk::generateProposalDensity(const VectorXd& x, const VectorXd& z, const MatrixXd& A, const VectorXd& b){
    generateHessian(x, A, b);
    VectorXd d = generateGaussianRV(x.rows());
    return sqrt(hess.determinant()) * exp(term_density * localNorm(x - z, hess));
}

void BarrierWalk::generateSample(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    generateHessian(x, A, b);
    MatrixXd matrix = hess.inverse().sqrt();
    VectorXd direction = generateGaussianRV(x.rows());
    z = x + term_sample * (matrix * direction);
}

void BarrierWalk::printType(){
    cout << "Generic Barrier" << endl;
}

MatrixXd BarrierWalk::generateCompleteWalk(const int num_steps, VectorXd& x, const MatrixXd& A, const VectorXd& b){
    MatrixXd results = MatrixXd::Zero(num_steps, A.cols());
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);
    float one = 1.0;

    float constant = (R * R)/b.rows();
    float td = (-0.5 / constant);
    float ts = sqrt(constant);
    BarrierWalk::setTs(ts);
    BarrierWalk::setTd(td);

    for(int i = 0; i < num_steps; i++){
        generateSample(x, A, b);
        if(inPolytope(z, A, b)){
            float g_x_z = generateProposalDensity(x, z, A, b);
            float g_z_x = generateProposalDensity(z, x, A, b);
            float alpha = min(one, g_z_x/g_x_z);
            float val = dis(gen);
            x = val < alpha ? z : x;
        }
        results.row(i) = x.transpose();
    }
    return results;
}