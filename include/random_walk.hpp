#ifndef RAND_H
#define RAND_H
#include "common.hpp"

class RandomWalk{
    public: 
        MatrixXd A{};
        VectorXd b{};
        float r {};

    RandomWalk(){}

    VectorXd vectPow(VectorXd& x, float alpha);
    virtual void initialize(MatrixXd A_p, VectorXd b_p, float r_p);
    bool acceptReject(VectorXd& z, MatrixXd& A, VectorXd& b);
    VectorXd generateGaussianRVNorm(int d);
    virtual MatrixXd generateCompleteWalk(int num_steps, VectorXd x);
};

#endif