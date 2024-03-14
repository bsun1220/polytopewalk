
#ifndef JOHNWALK_HPP
#define JOHNWALK_HPP

#include "BarrierWalk.hpp"

class JohnWalk: public BarrierWalk{

    public:
         /**
         * @brief constructor for John Walk class
         * @param r spread parameter
         * @param thin thin constant
         * @param g_lim gradient descent norm limit
         * @param step_size size of gradient descent step
         * @param max_iter maximum number of iterations in gradient descent
         */
        JohnWalk(double r, int thin = 1, double g_lim = 0.01, double step_size = 0.1, int max_iter = 100) : STEPSIZE(step_size), MAXITER(max_iter), GRADLIM(g_lim), BarrierWalk(r, thin){


        }

        /**
         * @brief print john walk
         * @return void
         */
        void printType() override;

        /**
         * @brief generates John weight by solving convex optimization problem
         * @param x point in polytope to generate DikinLS weight
         * @param A polytope matrix
         * @param b polytope matrix
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
        const double MAXITER;

        /**
         * @brief stops gradient descent if it reaches under this number
         */
        const double GRADLIM;

        /**
         * @brief set Dist Term for John Walk
         * @param d (dimension)
         * @param n (number of constraints)
         * @return void
         */
        void setDistTerm(int d, int n) override;

};

#endif