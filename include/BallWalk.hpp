#include "RandomWalk.hpp"

class BallWalk: public RandomWalk{
    
    public:
        BallWalk() : RandomWalk() {}

        /**
         * @brief Initialize values (because prior to Reducer, it is unknown what these values are)
         * @param A_p (Matrix for polytope) Ax <= b
         * @param b_p (Vector for polytope) Ax <= b
         * @param r_p general indicator of spread
         * @return void
         */
        void initialize(MatrixXd A_p, VectorXd b_p, float r_p){
            A = A_p;
            b = b_p;
            r = r_p;
        }

        /**
         * @brief Generate values from the walk
         * @param num_steps number of steps wanted to take
         * @param x initial starting point
         * @return Matrix
         */
        MatrixXd generateCompleteWalk(int num_steps, VectorXd& x);
        
        /**
         * @brief print general type 
         * @return void
         */
        void printType();


};