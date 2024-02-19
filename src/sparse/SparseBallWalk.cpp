#include "SparseBallWalk.hpp"

MatrixXd SparseBallWalk::generateCompleteWalk(
    const int num_steps, 
    const VectorXd& init, 
    const SparseMatrixXd& A, 
    const VectorXd& b, 
    int k
){
    MatrixXd results = MatrixXd::Zero(num_steps, A.cols());

    SparseLU<SparseMatrixXd> A_solver (A * A.transpose());
    SparseMatrixXd I = SparseMatrixXd(VectorXd::Ones(A.cols()).asDiagonal());

    VectorXd x = init;
    for (int i = 0; i < num_steps; i++){
        VectorXd rand = generateGaussianRV(A.cols()); 
        VectorXd z;
        z = A * rand; 
        z = rand - A.transpose() * A_solver.solve(z);
        z /= z.norm(); 
        z = R * z + x; 

        if (inPolytope(z, k)){
            x = z;
        }
        results.row(i) = x; 
    }
    return results; 
}