#ifndef CONSBARRIERWALK_HPP
#define CONSBARRIERWALK_HPP
#include "ConstraintWalk.hpp"

class ConstraintBarrierWalk : public ConstraintWalk{

    public:
        ConstraintBarrierWalk(const double err, const double r) : R(r), ConstraintWalk(err) {}

        virtual SparseMatrixXd generateG(
            const VectorXd& x, 
            const SparseMatrixXd& A,
            int k
        );

        MatrixXd generateCompleteWalk(
            const int num_steps, 
            const VectorXd& init, 
            const SparseMatrixXd& A, 
            const VectorXd& b, 
            int k
        ) override;
        
        virtual void setDistTerm(int d, int n); 
    
    protected:
        double DIST_TERM; 
        double R; 

        VectorXd generateSample(
            const VectorXd& x, 
            const SparseMatrixXd& A, 
            int k
        ); 
        
        double generateProposalDensity(
            const VectorXd& x, 
            const VectorXd& z, 
            const SparseMatrixXd& A, 
            int k
        );
};

#endif




