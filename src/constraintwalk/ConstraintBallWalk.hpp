#ifndef CONSBALLWALK_HPP
#define CONSBALLWALK_HPP

#include "ConstraintWalk.hpp"

class ConstraintBallWalk : public ConstraintWalk{
    public:
        ConstraintBallWalk(const double r) : R(r), ConstraintWalk(0.0){}

        MatrixXd generateCompleteWalk(
            const int num_steps, 
            const VectorXd& init, 
            const SparseMatrixXd& A, 
            const VectorXd& b, 
            int k) override;
    
    protected:
        const double R;
};
#endif 