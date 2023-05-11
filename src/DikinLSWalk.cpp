#include "DikinLSWalk.hpp"

void DikinLSWalk::gradientDescent(VectorXd& x, float adj, int sim, float gl){

    float q = 2 * (1 + log(A.rows()));
    float alpha = 1 - (2/q);

    VectorXd w_i = VectorXd::Ones(A.rows()); 
    generateSlack(x);
    MatrixXd slack_mat = slack.asDiagonal().toDenseMatrix();
    MatrixXd A_x = slack_mat.colPivHouseholderQr().solve(A);

    VectorXd term2 = (0.5 - 1/q) * VectorXd::Ones(A.rows()); 

    MatrixXd W(A.rows(), A.rows()); 
    VectorXd term1a (A.rows());
    VectorXd term1b (A.rows());
    VectorXd term1(A.rows());
    VectorXd gradient (A.rows()); 
    VectorXd proposal (A.rows()); 

    while(sim--){
        W = vectPow(w_i, alpha).asDiagonal().toDenseMatrix();
        term1a = alpha * vectPow(w_i, alpha - 1);
    
        term1b = (A_x * (A_x.transpose() * W * A_x).inverse() * A_x.transpose()).diagonal();

        term1 = term1a.cwiseProduct(term1b);
        
        gradient = term1 - term2;
        if(gradient.norm() < gl){
            break;
        }
        proposal = w_i + adj * gradient;
        if(proposal.minCoeff() < 0){
            break; 
        }
        w_i = proposal;
    }

    weights = w_i.asDiagonal().toDenseMatrix();
    
}

void DikinLSWalk::generateWeight(VectorXd& x){
    gradientDescent(x, step_size, max_iter, grad_lim);
}

void DikinLSWalk::printType(){
    cout << "DikinLSWalk" << endl;
}