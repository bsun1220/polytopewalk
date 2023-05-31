#ifndef POLY_H
#define POLY_H
#include "Common.hpp"
#include "CentralPointFinder.hpp"
#include "FacialReduction.hpp"
#include "BallWalk.hpp"
#include "DikinWalk.hpp"
#include "DikinLSWalk.hpp"
#include "JohnWalk.hpp"
#include "HitRunWalk.hpp"
#include "VaidyaWalk.hpp"

 /**
 * @brief does complete polytope walk process
 * @param A part of polytope initialization {x | Ax <= b}
 * @param b part of polytope initialization {x | Ax <= b}
 * @param walk pre-specified walk choice
 * @param reducer pre-specified reducer choice
 * @param initializer pre-specified initializer choice
 * @return Matrix
 */
MatrixXd fullWalkRun(MatrixXd A, VectorXd b, int num_sim, RandomWalk* walk, Reducer* reducer, Initializer* initializer){
    problem_result fr = reducer->reduce(A, b);
    int x_dim = A.cols();
    if (fr.reduced){
        MatrixXd reduced_A = fr.reduced_A;
        MatrixXd reduced_b = fr.reduced_b; 
        VectorXd pb = fr.b_tilde;
        MatrixXd M_inv = fr.M.inverse();

        VectorXd x = initializer->getInitialPoint(reduced_A, reduced_b);
        MatrixXd results = walk->generateCompleteWalk(num_sim, x, reduced_A, reduced_b);
        MatrixXd res(results.rows(), x_dim);
        for(int i = 0; i < results.rows(); i++){
            VectorXd val (results.cols() + pb.rows());
            VectorXd row = results.row(i);
            val << pb, row;
            res.row(i) = (M_inv * val).head(x_dim);
        }
        return res; 


    } else {

        VectorXd x = initializer->getInitialPoint(A, b);
        MatrixXd res = walk->generateCompleteWalk(num_sim, x, A, b);
        return res;
    }
}

#endif
