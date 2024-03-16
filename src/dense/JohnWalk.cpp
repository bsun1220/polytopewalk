#include "JohnWalk.hpp"

void JohnWalk::setDistTerm(int d, int n){
    DIST_TERM = R*R/(pow(d, 1.5));
}

void JohnWalk::generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b ){
    double alpha = 1 - 1/(log2(2 * A.rows() / A.cols()));
    double beta = (double)A.cols() / (2 * A.rows());

    VectorXd w_i = VectorXd::Ones(A.rows()); 
    generateSlack(x, A, b);
    DiagonalMatrix<double, Dynamic> slack_inv = slack.cwiseInverse().asDiagonal();

    MatrixXd A_x = slack_inv * A; 
    VectorXd term1 = VectorXd::Ones(A.rows());

    DiagonalMatrix<double, Dynamic> W;
    MatrixXd WAX (A.rows(), A.cols());
    VectorXd term2a (A.rows());
    VectorXd term2b (A.rows());
    VectorXd term2 (A.rows());
    VectorXd term3 (A.rows());
    VectorXd gradient (A.rows());

    VectorXd beta_ones = beta * VectorXd::Ones(A.rows());

    for(int i = 0; i < MAXITER; i++){
        W = VectorXd(w_i.array().pow(alpha * 0.5)).asDiagonal();
        WAX = W * A_x;
        term2a = w_i.cwiseInverse();
        term2b = (WAX * (WAX.transpose() * WAX).inverse()).cwiseProduct(WAX).rowwise().sum();

        term2 = term2a.cwiseProduct(term2b);
        term3 = beta * w_i.cwiseInverse();
        
        gradient = term1 - term2 - term3;
        if(gradient.norm() < GRADLIM){
            break;
        }
        w_i = (w_i - STEPSIZE * gradient).cwiseMax(beta_ones);
    }

    weights = w_i.asDiagonal();

}


void JohnWalk::printType(){
    cout << "John Walk" << endl;
}