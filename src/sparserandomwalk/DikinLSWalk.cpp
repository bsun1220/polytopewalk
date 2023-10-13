#include "DikinLSWalk.hpp"

void DikinLSWalk::generateWeight(const VectorXd& x, const SparseMatrixXd& A, const VectorXd& b){

    float q = 2 * (1 + log(A.rows()));
    float alpha = 1 - (2/q);

    VectorXd w_i = VectorXd::Ones(A.rows()); 
    generateSlack(x, A, b);
    SparseMatrixXd slack_mat = SparseMatrixXd(slack.asDiagonal());

    SimplicialCholesky<SparseMatrixXd> chol (slack_mat);
    SparseMatrixXd A_x = chol.solve(A);
    SparseMatrixXd A_x_T = A_x.transpose();

    VectorXd term2 = (0.5 - 1/q) * VectorXd::Ones(A.rows()); 

    SparseMatrixXd W(A.rows(), A.rows()); 
    VectorXd term1a (A.rows());
    VectorXd term1b (A.rows());
    VectorXd term1(A.rows());
    VectorXd gradient (A.rows()); 
    VectorXd proposal (A.rows()); 

    for(int i = 0; i < MAXITER; i++){
        W = SparseMatrixXd(vectPow(w_i, alpha).asDiagonal());
        term1a = alpha * vectPow(w_i, alpha - 1);
    
        SimplicialCholesky<SparseMatrixXd> chol2 (A_x_T * W * A_x);
        SparseMatrixXd res = A_x * chol2.solve(A_x_T);
        term1b = res.diagonal();

        term1 = term1a.cwiseProduct(term1b);
        
        gradient = term1 - term2;
        if(gradient.norm() < GRADLIM){
            break;
        }
        proposal = w_i + STEPSIZE * gradient;
        if(proposal.minCoeff() < 0){
            break; 
        }
        w_i = proposal;
    }

    weights = SparseMatrixXd(w_i.asDiagonal());
    
}

void DikinLSWalk::printType(){
    cout << "DikinLSWalk" << endl;
}