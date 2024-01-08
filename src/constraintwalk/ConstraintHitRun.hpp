#ifndef CONSHITRUN_HPP
#define CONSHITRUN_HPP

#include "ConstraintWalk.hpp"

class ConstraintHitAndRun : public ConstraintWalk{
    public:
        ConstraintHitAndRun(const double err, const double r) : R(r), ConstraintWalk(err) {}

        MatrixXd generateCompleteWalk(
            const int num_steps, 
            const VectorXd& init, 
            const SparseMatrixXd& A, 
            const VectorXd& b, 
            int k) override;
        
        double binarySearch(
            VectorXd direction, 
            VectorXd& x,
            int k);
    
    protected:
        const double R;
};

#endif