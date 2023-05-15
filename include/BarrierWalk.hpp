#ifndef BARRIER_H
#define BARRIER_H

#include "RandomWalk.hpp"

class BarrierWalk : public RandomWalk{
    public:
        
        BarrierWalk(){
            term_sample = 0;
            term_density = 0;
        }

        /**
         * @brief Initialize values (because prior to Reducer, it is unknown what these values are)
         * @param A_p (Matrix for polytope) Ax <= b
         * @param b_p (Vector for polytope) Ax <= b
         * @param r_p general indicator of spread
         * @return void
         */
        virtual void initialize(MatrixXd A_p, VectorXd b_p, float r);

        /**
         * @brief Generate values from the walk
         * @param num_steps number of steps wanted to take
         * @param x initial starting point
         * @return Matrix
         */
        MatrixXd generateCompleteWalk(const int num_steps, VectorXd& x);
    
    protected:
        /**
         * @brief term_sample for coefficient in sample function
         */
        float term_sample{};

        /**
         * @brief term_sample for coefficient in density function
         */
        float term_density{};

        /**
         * @brief represents global variable b - Ax
         */
        VectorXd slack{}; 

        /**
         * @brief Hessian Matrix from global variable from generateHessian
         */
        MatrixXd hess{};

        /**
         * @brief new proposal point generated from generateSample function
         */
        VectorXd z{};

        /**
         * @brief weights generated from generateWeights function
         */
        MatrixXd weights{};

        /**
         * @brief set term sample by indicated value
         * @param a
         * @return void
         */
        void setTs(float a);

        /**
         * @brief set term density by indicated value
         * @param b
         * @return void
         */
        void setTd(float b);

        /**
         * @brief generates a gaussian random vector with d dimension
         * @param d dimension
         * @return Vector
         */
        VectorXd generateGaussianRV(int d);

        /**
         * @brief generates b - Ax (called slack) and 
         * makes global variable slack equal to it
         * @param x
         * @return void
         */
        void generateSlack(VectorXd& x);

        /**
         * @brief calculates distance weighted by Hessian matrix m
         * @param m Weighted Hessian Matrix
         * @param v vector to be measured
         * @return float 
         */
        float localNorm(VectorXd v, MatrixXd& m);

        /**
         * @brief generate weights when calculating Hessian matrix
         * @param x point in polytope to generate weight
         * @return void (update global variable weights)
         */
        virtual void generateWeight(VectorXd& x);

        /**
         * @brief generate weights when calculating Hessian matrix
         * @param x point in polytope to generate weight
         * @return void (update global variable weights)
         */
        virtual void printType();

        /**
         * @brief print general type 
         * @return void
         */
        void generateHessian(VectorXd& x);

        /**
         * @brief generate proposal density around as a Multivariate Gaussian N(x, f(Hessian(x)))
         * where function f: R^N -> R^NxN varies depending on walk type for vector z
         * @param x centered point around ellipsoid
         * @param z new proposal point
         * @return float (representing density)
         */
        float generateProposalDensity(VectorXd& x, VectorXd& z);

        /**
         * @brief generates a point drawn from a Multivariate Gaussian N(x, f(Hessian(x)))
         * @param x centered point in the polytope
         * @return void (updates global variable z)
         */
        void generateSample(VectorXd& x);
};

#endif