
#include "utils/FacialReduction.hpp"

VectorXd generateGaussianRV(int d){
    VectorXd v(d);
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dis(0.0, 1.0);
    for(int i = 0; i < d; i++){
        v(i) = dis(gen);
    }
    return v;
}

int main(){

    MatrixXd A (1,3);
    A << 1, 1, 1;
    VectorXd b (1);
    b << 1;

    HouseholderQR <MatrixXd> qr(A.cols(), A.rows());
    qr.compute(A.transpose());
    MatrixXd Q = qr.householderQ();
    MatrixXd R  = qr.matrixQR().triangularView<Eigen::Upper>();
    int d = R.rows();
    int n = R.cols();

    MatrixXd newR = R.block(0, 0, R.cols(), R.cols());

    VectorXd z1 = newR.transpose().inverse() * b;
    
    MatrixXd Q1 = Q.block(0, 0, Q.rows(), n);
    MatrixXd Q2 = Q.block(0, n, Q.rows(), d - n);
    
    MatrixXd reduced_A = -1 * Q2;
    VectorXd reduced_b = (Q1 * z1);

    VectorXd l (3); 
    l << 0, 0, 0;

    VectorXd r (3);
    r << 1, 1, 1;

    VectorXd x (3);
    x << 0.33, 0.33, 0.34;

    MatrixXd I = MatrixXd::Identity(3, 3);
    MatrixXd Ap (6, 3);
    Ap << -I, I;
    VectorXd bp (6);
    bp << -l, r;

    MatrixXd slack_inv = (bp - (Ap * x)).cwiseInverse().asDiagonal().toDenseMatrix();
    MatrixXd G = Ap.transpose() * slack_inv * slack_inv * Ap;

    MatrixXd G_inv = G.inverse();
    MatrixXd G_inv_sqrt = G.inverse().sqrt();

    MatrixXd P = G_inv_sqrt * A.transpose() * (A * G_inv * A.transpose()).inverse() * A * G_inv_sqrt;

    MatrixXd M_inv = G_inv_sqrt * (I - P) * G_inv_sqrt;

    VectorXd rand = generateGaussianRV(3);

    VectorXd v = x + M_inv * rand;
    cout << v.sum() << endl;
    cout << v << endl;



}
