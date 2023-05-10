#ifndef BARRIER_H
#define BARRIER_H

#include "random_walk.hpp"

class BarrierWalk : public RandomWalk{
    public:
        MatrixXd A{};
        VectorXd b{};
        float term_sample{};
        float term_density{};

        VectorXd slack{}; 
        MatrixXd hess{};
        VectorXd z{};
        MatrixXd weights{};
        
        BarrierWalk(){
            term_sample = 0;
            term_density = 0;
        }

        virtual void initialize(MatrixXd A_p, VectorXd b_p, float r);

        void setTs(float a);
        void setTd(float b);

        VectorXd generateGaussianRV(int d);
        void generateSlack(VectorXd& x);
        float localNorm(VectorXd v, MatrixXd& m);

        virtual void generateWeight(VectorXd& x);
        virtual void printType();

        void generateHessian(VectorXd& x);
        float generateProposalDensity(VectorXd& x, VectorXd& z);
        void generateSample(VectorXd& x);

        MatrixXd generateCompleteWalk(int num_steps, VectorXd& x);
};



#endif