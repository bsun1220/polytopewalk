#include "BarrierWalk.hpp"

void BarrierWalk::setDistTerm(int d, int n){
    DIST_TERM = R*R/n;
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

double BarrierWalk::localNorm(VectorXd v, const MatrixXd& m){
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

double BarrierWalk::generateProposalDensity(const VectorXd& x, const VectorXd& z, const MatrixXd& A, const VectorXd& b){
    generateHessian(x, A, b);
    VectorXd d = generateGaussianRV(x.rows());

    LLT<MatrixXd> cholesky(hess);
    MatrixXd L = cholesky.matrixL();
    double det = L.diagonal().array().log().sum(); 
    double dist = -(0.5/DIST_TERM) * localNorm(x - z, hess);
    return det + dist; 
}

void BarrierWalk::generateSample(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    generateHessian(x, A, b);
    LLT<MatrixXd> cholesky(hess);
    MatrixXd L = cholesky.matrixL();
    FullPivLU<MatrixXd> lu(L);
    VectorXd direction = generateGaussianRV(x.rows());
    z = x + sqrt(DIST_TERM) * (lu.solve(direction));
}

MatrixXd BarrierWalk::generateCompleteWalk(const int num_steps, VectorXd& x, const MatrixXd& A, const VectorXd& b){
    MatrixXd results = MatrixXd::Zero(num_steps, A.cols());
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);
    double one = 1.0;

    setDistTerm(A.cols(), A.rows());

    for(int i = 0; i < num_steps; i++){
        generateSample(x, A, b);
        if(inPolytope(z, A, b)){
            double g_x_z = generateProposalDensity(x, z, A, b);
            double g_z_x = generateProposalDensity(z, x, A, b);
            double alpha = min(one, exp(g_z_x-g_x_z));
            double val = dis(gen);
            x = val < alpha ? z : x;
        }
        results.row(i) = x.transpose();
    }
    return results;
}