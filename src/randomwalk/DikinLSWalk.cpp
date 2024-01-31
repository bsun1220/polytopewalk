#include "DikinLSWalk.hpp"

void DikinLSWalk::setDistTerm(int d, int n){
    double q = 2.0 * (1.0 + log(n));
    double term = (1.0 + q) * (1.0 + q * q);
    DIST_TERM = R*R/term; 
}

void DikinLSWalk::generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b){

    float q = 2 * (1 + log(A.rows()));
    float alpha = 1 - (2/q);

    VectorXd w_i = VectorXd::Ones(A.rows()); 
    generateSlack(x, A, b);
    MatrixXd slack_mat = slack.asDiagonal().toDenseMatrix();
    MatrixXd A_x = slack_mat.colPivHouseholderQr().solve(A);

    VectorXd term1 = (0.5 - 1/q) * VectorXd::Ones(A.rows()); 

    MatrixXd W(A.rows(), A.rows()); 
    VectorXd term2a (A.rows());
    VectorXd term2b (A.rows());
    VectorXd term2(A.rows());
    VectorXd gradient (A.rows()); 
    VectorXd proposal (A.rows()); 
    VectorXd term3 (A.rows());
    double beta = = (double)A.cols() / (2 * A.rows());

    for(int i = 0; i < MAXITER; i++){
        W = vectPow(w_i, alpha).asDiagonal().toDenseMatrix();
        term2a = alpha * vectPow(w_i, alpha - 1);
    
        term2b = (A_x * (A_x.transpose() * W * A_x).inverse() * A_x.transpose()).diagonal();

        term2 = term2a.cwiseProduct(term2b);
        term3 = beta * w_i.cwiseInverse();
        
        gradient = term1 - term2 - term3; 
        if(gradient.norm() < GRADLIM){
            break;
        }
        w_i = w_i - STEPSIZE * gradient;
    }

    weights = w_i.asDiagonal().toDenseMatrix();
    
}

void DikinLSWalk::printType(){
    cout << "DikinLSWalk" << endl;
}