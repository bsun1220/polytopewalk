#include "Common.hpp"
#include "CentralPointFinder.hpp"
#include "FacialReduction.hpp"
#include "BallWalk.hpp"
#include "DikinWalk.hpp"
#include "DikinLSWalk.hpp"
#include "JohnWalk.hpp"
#include "HitRunWalk.hpp"
#include "VaidyaWalk.hpp"

VectorXd convertBack(VectorXd z, VectorXd& pb, MatrixXd& M_inv, int x_dim){
    VectorXd val (z.rows() + pb.rows());
    val << pb, z;
    return (M_inv * val).head(x_dim);
};

MatrixXd fullWalkRun(MatrixXd A, VectorXd b, float r, int num_sim, RandomWalk& walk, Reducer& reducer, Initializer& initializer){
    problem_result fr = reducer.reduce(A, b);
    int x_dim = A.cols();
    if (fr.reduced){
        MatrixXd reduced_A = fr.reduced_A;
        MatrixXd reduced_b = fr.reduced_b; 
        VectorXd pb = fr.b_tilde;
        MatrixXd M_inv = fr.M.inverse();

        walk.initialize(reduced_A, reduced_b, r);
        VectorXd x = initializer.getInitialPoint(reduced_A, reduced_b);
        MatrixXd results = walk.generateCompleteWalk(num_sim, x);

        MatrixXd res(results.rows(), x_dim);
        for(int i = 0; i < results.rows(); i++){
            res.row(i) = convertBack(results.row(i), pb, M_inv, x_dim);
        }
        return res; 


    } else {
        walk.initialize(A, b, r);

        VectorXd x = initializer.getInitialPoint(A, b);
        MatrixXd res = walk.generateCompleteWalk(num_sim, x);
        return res;
    }
}
