#ifndef CONSVAIDYA_HPP
#define CONSVAIDYA_HPP

#include "ConstraintBarrierWalk.hpp"

class ConstraintVaidyaWalk : public ConstraintBarrierWalk{

    public:
        ConstraintVaidyaWalk(const double err, const double r) : ConstraintBarrierWalk(err, r) {}
    
    protected:
        SparseMatrixXd generateG(
            const VectorXd& x, 
            const SparseMatrixXd& A,
            int k
        ) override; 

        void setDistTerm(int d, int n) override;

};

#endif