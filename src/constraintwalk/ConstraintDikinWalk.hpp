#ifndef CONSDIKIN_HPP
#define CONSDIKIN_HPP

#include "ConstraintBarrierWalk.hpp"

class ConstraintDikinWalk : public ConstraintBarrierWalk{

    public:
        ConstraintDikinWalk(const double err, const double r) : ConstraintBarrierWalk(err, r) {}
    
    protected:
        SparseMatrixXd generateG(
            const VectorXd& x, 
            const SparseMatrixXd& A,
            int k
        ) override; 

        void setDistTerm(int d, int n) override;

};

#endif

