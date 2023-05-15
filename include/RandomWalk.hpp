#ifndef RAND_H
#define RAND_H
#include "Common.hpp"

class RandomWalk{

    public:
    
        RandomWalk(){}

        /**
         * @brief Initialize values (because prior to Reducer, it is unknown what these values are)
         * @param A_p (Matrix for polytope) Ax <= b
         * @param b_p (Vector for polytope) Ax <= b
         * @param r_p general indicator of spread
         * @return void
         */
        virtual void initialize(MatrixXd A_p, VectorXd b_p, float r_p);

        /**
         * @brief Generate values from the walk
         * @param num_steps number of steps wanted to take
         * @param x initial starting point
         * @return Matrix
         */
        virtual MatrixXd generateCompleteWalk(const int num_steps, VectorXd& x);

    protected: 

        /**
         * @brief matrix for x {x | Ax <= b}
         */
        MatrixXd A;

        /**
         * @brief vector for x {x | Ax <= b}
         */
        VectorXd b;

        /**
         * @brief general spread hyperparameter
         */
        float r;

        /**
         * @brief return elementwise [x_1^alpha,...,a_n^alpha]
         * @param x
         * @param alpha
         * @return vector
         */
        VectorXd vectPow(VectorXd& x, const float alpha);

        /**
         * @brief checks Az <= b
         * @param z
         * @param A
         * @param x
         * @return void
         */
        bool acceptReject(VectorXd& z, MatrixXd& A, VectorXd& b);

        /**
         * @brief returns normalized Gaussian vector of dimension d
         * @param d
         * @return vector 
         */
        VectorXd generateGaussianRVNorm(const int d);

        /**
         * @brief prints unique identifier of the walk
         * @return void
         */
        virtual void printType();
};

#endif