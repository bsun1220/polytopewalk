#include <randomwalk/PolytopeWalk.hpp>


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

MatrixXd generateSlack(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    MatrixXd slack = (b - (A * x));
    return slack;
}

MatrixXd generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    int d = b.rows();
    MatrixXd weights = MatrixXd::Identity(d, d);
    return weights;
}

MatrixXd generateHessian(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    MatrixXd weights = generateWeight(x, A, b);
    MatrixXd slack = generateSlack(x, A, b);
    MatrixXd slack_inv = slack.cwiseInverse().asDiagonal().toDenseMatrix();
    return A.transpose() * slack_inv * weights * slack_inv * A;

}

float localNorm(VectorXd v, const MatrixXd& m){
    return ((v.transpose() * m) * v)(0);
}

float generateProposalDensity(const VectorXd& x, const VectorXd& z, const MatrixXd& A, const VectorXd& b){
    MatrixXd hess = generateHessian(x, A, b);
    VectorXd d = generateGaussianRV(x.rows());
    return sqrt(hess.determinant()) * exp(localNorm(x - z, hess));
}

VectorXd generateSample(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    MatrixXd hess = generateHessian(x, A, b);
    MatrixXd matrix = hess.inverse().sqrt();
    VectorXd direction = generateGaussianRV(x.rows());
    VectorXd z = x + (matrix * direction);
    return z;
}

int main(){

    MatrixXd A (4,2);
    A << 1, 0, -1, 0, 0, 1, 0, -1;
    VectorXd b(4);
    b << 1, 0, 1, 0;
    VectorXd x(2);
    x << 0.5, 0.5;
    VectorXd z(2);
    z << 0.75, 0.75;

    /*
    MatrixXd A (4,2);
    A << -1, 0, 0, -1, 1, 1, -1, 1;
    VectorXd b(4);
    b << 0, 0, 1, 0;
    VectorXd x(2);
    x << 0.25, 0.5;
    VectorXd z(2);
    z << 0.4, 0.3;*/


    MatrixXd hess = generateHessian(x, A, b);
    SimplicialLLT<SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> cholesky;
    cout << hess << endl;
    SparseMatrix<double> shess = hess.sparseView();

    cholesky.analyzePattern(shess);
    cholesky.factorize(shess);

    /**
    std::cout << Eigen::MatrixXd(shess) << std::endl;
    std::cout << Eigen::MatrixXd(cholesky.matrixU()) << std::endl;
    std::cout << Eigen::MatrixXd(cholesky.matrixL()) << std::endl;
    **/
   
    MatrixXd L = MatrixXd(cholesky.matrixL());
    cout << "----" << endl;
    cout << pow(L.diagonal().prod(),2)  << endl;
    cout << hess.determinant() << endl;
    cout << "----" << endl;
    MatrixXd I = Eigen::MatrixXd::Identity(2,2);
    
    cout << cholesky.solve(I) << endl;

    cout << hess.inverse() << endl;
    cout << "----" << endl;
    cout << hess.inverse().sqrt() << endl;
    cout << cholesky.solve(L) << endl;

}