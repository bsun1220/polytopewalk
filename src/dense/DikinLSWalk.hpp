
#ifndef DIKINLSWALK_HPP
#define DIKINLSWALK_HPP

#include "BarrierWalk.hpp"

class DikinLSWalk: public BarrierWalk{

    public:
        /**
         * @brief constructor for weighted Dikin Walk class
         * @param r spread parameter
         * @param g_lim gradient descent norm limit
         * @param step_size size of gradient descent step
         * @param max_iter maximum number of iterations in gradient descent
         */
        DikinLSWalk(const double r, const double g_lim, const double step_size, const int max_iter) : STEPSIZE(step_size), MAXITER(max_iter), GRADLIM(g_lim), BarrierWalk(r){
            
        }

        /**
         * @brief print dikinls
         * @return void
         */
        void printType() override;

        /**
         * @brief generate weights when calculating Hessian matrix
         * @param x point in polytope to generate DikinLS weight
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         * @return void (update global variable weights)
         */
        void generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b) override;
    
    protected:
        /**
         * @brief step size for gradient descent
         */
        const double STEPSIZE;

        /**
         * @brief max number of iterations in gradient descent
         */
        const int MAXITER;

        /**
         * @brief stops gradient descent if it reaches under this number
         */
        const double GRADLIM;

        /**
         * @brief set Distribution Constant
         * @param d (dimension)
         * @param n (number of constraints)
         * @return void
         */
        void setDistTerm(int d, int n) override;
};

#endif