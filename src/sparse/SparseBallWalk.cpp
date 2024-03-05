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
    int d = A.cols() - A.rows();
    int total = num_steps * THIN;
    for (int i = 1; i <= total; i++){
        VectorXd rand = generateGaussianRV(A.cols()); 
        VectorXd z;
        z = A * rand; 
        z = rand - A.transpose() * A_solver.solve(z);
        z /= z.norm(); 
        z = R/sqrt(d) * z + x; 

        if (inPolytope(z, k)){
            x = z;
        }
        if (i % THIN == 0){
            results.row((int)i/THIN - 1) = x; 
        }
    }
    return results; 
}