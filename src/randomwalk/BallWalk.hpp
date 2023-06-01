
#ifndef BALLWALK_HPP
#define BALLWALK_HPP

#include "RandomWalk.hpp"

class BallWalk: public RandomWalk{
    

    public:

        /**
         * @brief Ball Walk Implementation
         * @param r_p radius of ball hyper parameterx
         */
        BallWalk(const float r_p) : R(r_p), RandomWalk() {
            
        }

        /**
         * @brief Generate values from the walk
         * @param num_steps number of steps wanted to take
         * @param x initial starting point
         * @param A polytope
         * @param b polytope
         * @return Matrix
         */
        MatrixXd generateCompleteWalk(const int num_steps, VectorXd& x, const MatrixXd& A, const VectorXd& b) override;
        
        /**
         * @brief print general type 
         * @return void
         */
        void printType() override;
    
    protected:
        const float R;


};

#endif