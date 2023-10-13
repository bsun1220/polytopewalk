#include "JohnWalk.hpp"
#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>

void JohnWalk::generateWeight(const VectorXd& x, const SparseMatrixXd& A, const VectorXd& b ){
    float alpha = 1 - 1/(log2(2 * A.rows() / A.cols()));
    float beta = (double)A.cols() / (2 * A.rows());

    VectorXd w_i = VectorXd::Ones(A.rows()); 
    generateSlack(x, A, b);

    SparseMatrixXd slack_mat = SparseMatrixXd(slack.asDiagonal());

    SimplicialCholesky<SparseMatrixXd> chol (slack_mat);
    SparseMatrixXd A_x = chol.solve(A);
    SparseMatrixXd A_x_T = A_x.transpose();

    VectorXd term1 = VectorXd::Ones(A.rows());
    
    SparseMatrixXd W(A.rows(), A.rows());
    VectorXd term2a (A.rows());
    VectorXd term2b (A.rows());
    VectorXd term2 (A.rows());
    VectorXd term3 (A.rows());
    VectorXd gradient (A.rows());

    for(int i = 0; i < MAXITER; i++){
        W = SparseMatrixXd(vectPow(w_i, alpha).asDiagonal());

        term2a = vectPow(w_i, alpha - 1);

        SimplicialCholesky<SparseMatrixXd> chol2 (A_x_T * W * A_x);
        SparseMatrixXd res = A_x * chol2.solve(A_x_T);
        term2b = res.diagonal();

        term2 = term2a.cwiseProduct(term2b);
        term3 = beta * w_i.cwiseInverse();
        
        gradient = term1 - term2 - term3;
        if(gradient.norm() < GRADLIM){
            break;
        }
        w_i = w_i - STEPSIZE * gradient;
    }

    weights = SparseMatrixXd(w_i.asDiagonal());

}


void JohnWalk::printType(){
    cout << "John Walk" << endl;
}