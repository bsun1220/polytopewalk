#ifndef CONSTRAINTWALK_HPP
#define CONSTRAINTWALK_HPP
#include "Common.hpp"

class ConstraintWalk{

    public:
        ConstraintWalk(const double err) : ERR(err){}
    
        virtual MatrixXd generateCompleteWalk(
            const int num_steps, 
            const VectorXd& init, 
            const SparseMatrixXd& A, 
            const VectorXd& b, 
            int k);
        
    protected:
        bool inPolytope(const VectorXd& z, int k);
        bool inApproxPolytope(const VectorXd& z, int k, double err);
        VectorXd generateGaussianRV(const int d);
        const double ERR; 
};

#endif