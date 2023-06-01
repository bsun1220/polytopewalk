
#ifndef HITRUN_HPP
#define HITRUN_HPP

#include "RandomWalk.hpp"

class HitAndRunWalk: public RandomWalk{

    public:
        /**
         * @brief Hit and Run implementation constructor
         * @param err_p error hyperparameter
         * @param r spread hyperparamter
         */
        HitAndRunWalk(const float err_p, const float r) : ERR(err_p), R(r), RandomWalk() {

        }

        /**
         * @brief Generate values from the walk
         * @param num_steps number of steps wanted to take
         * @param x initial starting point
         * @param A polytope matrix
         * @param b polytope matrix
         * @return Matrix
         */
        MatrixXd generateCompleteWalk(const int num_steps, VectorXd& x, const MatrixXd& A, const VectorXd& b) override;

         /**
         * @brief print general type 
         * @return void
         */
        void printType() override;
    
    protected:
        /**
         * @brief relative error of the binary search operation
         */
        const float ERR;

        /**
         * @brief initial starting value
         */
        const float R;

        /**
         * @brief get distance between vectors x and y
         * @param x
         * @param y
         * @return double 
         */
        double distance(VectorXd& x, VectorXd&y);

        /**
         * @brief runs binary search to find a suitable chord intersection with the polytope
         * @param direction (random direction variable)
         * @param x (starting point)
         * @param A polytope matrix
         * @param b polytope vector
         * @return double 
         */
        double binarySearch(VectorXd direction, VectorXd& x, const MatrixXd& A, const VectorXd& b);

};

#endif