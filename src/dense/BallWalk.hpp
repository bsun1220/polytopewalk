
#ifndef BALLWALK_HPP
#define BALLWALK_HPP

#include "RandomWalk.hpp"

class BallWalk: public RandomWalk{
    

    public:

        BallWalk(const double r) : R(r), RandomWalk() {
            
        }

        /**
         * @brief Generate values from Ball Walk
         * @param num_steps number of steps wanted to take
         * @param x initial starting point
         * @param A polytope matrixd (Ax <= b)
         * @param b polytope vector (Ax <= b)
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
         * @brief spread parameter
         */
        const double R;


};

#endif