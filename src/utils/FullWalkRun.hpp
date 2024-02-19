#include "FacialReduction.hpp"
#include "DenseCenter.hpp"
#include "SparseCenter.hpp"
#include "dense/DikinWalk.hpp"
#include "dense/DikinLSWalk.hpp"
#include "dense/JohnWalk.hpp"
#include "dense/VaidyaWalk.hpp"
#include "dense/HitRun.hpp"
#include "dense/BallWalk.hpp"
#include "sparse/SparseDikinWalk.hpp"
#include "sparse/SparseDikinLSWalk.hpp"
#include "sparse/SparseJohnWalk.hpp"
#include "sparse/SparseVaidyaWalk.hpp"
#include "sparse/SparseBallWalk.hpp"
#include "sparse/SparseHitRun.hpp"

MatrixXd denseFullWalkRun(SparseMatrixXd A, VectorXd b, int k, int num_sim, RandomWalk* walk){
    FacialReduction fr;
    res fr_result = fr.reduce(A, b, k, false);
    DenseCenter init; 
    VectorXd x = init.getInitialPoint(fr_result.dense_A, fr_result.dense_b);
    MatrixXd steps = walk->generateCompleteWalk(num_sim, x, fr_result.dense_A, fr_result.dense_b);
    MatrixXd res(num_sim, A.cols());
    for(int i = 0; i < num_sim; i++){
        VectorXd val (steps.cols() + fr_result.z1.rows());
        VectorXd row = steps.row(i);
        val << fr_result.z1, row;
        res.row(i) = (fr_result.Q * val).head(A.cols());
    }
    return res; 
}

MatrixXd sparseFullWalkRun(SparseMatrixXd A, VectorXd b, int k, int num_sim, SparseRandomWalk* walk){
    FacialReduction fr;
    res fr_result = fr.reduce(A, b, k, true);
    int new_k = fr_result.sparse_A.rows() - (A.rows() - k);
    SparseCenter init; 
    VectorXd x = init.getInitialPoint(fr_result.sparse_A, fr_result.sparse_b, new_k);
    MatrixXd steps = walk->generateCompleteWalk(num_sim, x, fr_result.sparse_A, fr_result.sparse_b, new_k);
    MatrixXd res(num_sim, A.cols());
    for(int i = 0; i < num_sim; i++){
        res.row(i) = fr_result.saved_V * steps.row(i).transpose();
    }
    return res; 
}