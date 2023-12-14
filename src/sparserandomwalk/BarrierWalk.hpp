#ifndef BARRIER_HPP
#define BARRIER_HPP

#include "RandomWalk.hpp"
#include "utils/LeverageScore.hpp"

class BarrierWalk : public RandomWalk{
    public:
        
        /**
         * @brief BarrierWalk class
         * @param rp spread parameter
         */
        BarrierWalk(const float rp) : R(rp), RandomWalk(){

        }

        /**
         * @brief generate weights when calculating Hessian matrix
         * @param x point in polytope to generate weight
         * @param A polytope matrix
         * @param b polytope vector
         * @return void (update global variable weights)
         */
        virtual void generateWeight(const VectorXd& x, const SparseMatrixXd& A, const VectorXd& b);

        /**
         * @brief Generate values from the walk
         * @param num_steps number of steps wanted to take
         * @param x initial starting point
         * @param A polytope matrix
         * @param b polytope vector
         * @return Matrix
         */
        MatrixXd generateCompleteWalk(const int num_steps, VectorXd& x, const SparseMatrixXd& A, const VectorXd& b) override;

        /**
         * @brief set Dist Term
         * @param d
         * @param n
         * @return void
         */
        virtual void setDistTerm(int d, int n);
    
    protected:

        /**
         * @brief spread parameter
         */
        const float R;

        /**
         * @brief normalized spread parameter
         */
        float DIST_TERM;

        /**
         * @brief represents global variable b - Ax
         */
        VectorXd slack{}; 

        /**
         * @brief Hessian Matrix from global variable from generateHessian
         */
        SparseMatrixXd hess{};

        /**
         * @brief new proposal point generated from generateSample function
         */
        VectorXd z{};

        /**
         * @brief weights generated from generateWeights function
         */
        SparseMatrixXd weights{};

        /**
         * @brief generates a gaussian random vector with d dimension
         * @param d dimension
         * @return Vector
         */
        VectorXd generateGaussianRV(int d);

        /**
         * @brief generates b - Ax (called slack) and 
         * makes global variable slack equal to it
         * @param x point
         * @param A polytope matrix
         * @param b polytope vector
         * @return void
         */
        void generateSlack(const VectorXd& x, const SparseMatrixXd& A, const VectorXd& b);

        /**
         * @brief calculates distance weighted by Hessian matrix m
         * @param m Weighted Hessian Matrix
         * @param v vector to be measured
         * @return float 
         */
        float localNorm(VectorXd v, const SparseMatrixXd& m);

        /**
         * @brief print general type 
         * @return void
         */
        void generateHessian(const VectorXd& x, const SparseMatrixXd& A, const VectorXd& b);

        /**
         * @brief generate proposal density around as a Multivariate Gaussian N(x, f(Hessian(x)))
         * where function f: R^N -> R^NxN varies depending on walk type for vector z
         * @param x centered point around ellipsoid
         * @param z new proposal point
        * @param A polytope matrix
         * @param b polytope vector
         * @return float (representing density)
         */
        float generateProposalDensity(const VectorXd& x, const VectorXd& z, const SparseMatrixXd&A, const VectorXd& b);

        /**
         * @brief generates a point drawn from a Multivariate Gaussian N(x, f(Hessian(x)))
         * @param x centered point in the polytope
         * @param A polytope matrix
         * @param b polytope vector
         * @return void (updates global variable z)
         */
        void generateSample(const VectorXd& x, const SparseMatrixXd& A, const VectorXd& b);
};

#endif