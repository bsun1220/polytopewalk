
#ifndef JOHNWALK_HPP
#define JOHNWALK_HPP

#include "BarrierWalk.hpp"

class JohnWalk: public BarrierWalk{

    public:
        /**
         * @brief constructor for weighted John walk class
         * @param ss 
         * @param mi 
         * @param gl 
         * @param r
         */
        JohnWalk(const float ss, const float mi, const float gl, const float r) : STEPSIZE(ss), MAXITER(mi), GRADLIM(gl), BarrierWalk(r){


        }

        /**
         * @brief print john walk
         * @return void
         */
        void printType() override;
    
    protected:
        /**
         * @brief step size for gradient descent
         */
        const float STEPSIZE;

        /**
         * @brief max number of iterations in gradient descent
         */
        const float MAXITER;

        /**
         * @brief stops gradient descent if it reaches under this number
         */
        const float GRADLIM;

        /**
         * @brief generate weights when calculating Hessian matrix
         * @param x point in polytope to generate DikinLS weight
         * @param A polytope matrix
         * @param b polytope matrix
         * @return void (update global variable weights)
         */
        void generateWeight(const VectorXd& x, const SparseMatrixXd& A, const VectorXd& b) override;

};

#endif