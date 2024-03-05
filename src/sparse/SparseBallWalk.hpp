#ifndef CONSBALLWALK_HPP
#define CONSBALLWALK_HPP

#include "SparseRandomWalk.hpp"

class SparseBallWalk : public SparseRandomWalk{
    public:
        /**
         * @brief SparseBall Walk class
         * @param r spread parameter
         * @param thin thin parameter
         */
        SparseBallWalk(double r, int thin = 1) : R(r), SparseRandomWalk(0.0, thin){}

         /**
         * @brief Generate values from the Ball walk
         * @param num_steps number of steps wanted to take
         * @param init initial starting point
         * @param A polytope matrix 
         * @param b polytope vector
         * @param k k values >= 0 constraint
         * @return Matrix
         */
        MatrixXd generateCompleteWalk(
            const int num_steps, 
            const VectorXd& init, 
            const SparseMatrixXd& A, 
            const VectorXd& b, 
            int k) override;
    
    protected:
        /**
         * @brief spread parameter
         */
        const double R;
};
#endif 