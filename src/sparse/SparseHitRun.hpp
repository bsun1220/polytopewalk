#ifndef CONSHITRUN_HPP
#define CONSHITRUN_HPP

#include "SparseRandomWalk.hpp"

class SparseHitAndRun : public SparseRandomWalk{
    public:
        /**
         * @brief constructor for HitAndRun
         * @param err error constant
         * @param r spread parameter
         */
        SparseHitAndRun(const double err, const double r) : R(r), SparseRandomWalk(err) {}

         /**
         * @brief Generate values from the Hit and Run
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

        /**
         * @brief runs binary search to find a suitable chord intersection with the polytope
         * @param direction (random direction variable)
         * @param x (starting point)
         * @param k k values >= 0 constraint
         * @return double 
         */
        double binarySearch(
            VectorXd direction, 
            VectorXd& x,
            int k);
};

#endif