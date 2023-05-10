#include "random_walk.hpp"

class BallWalk: public RandomWalk{
    
    public:
        BallWalk() : RandomWalk() {}

        void initialize(MatrixXd A_p, VectorXd b_p, float r_p){
            A = A_p;
            b = b_p;
            r = r_p;
        }

        MatrixXd generateCompleteWalk(int num_steps, VectorXd x);

};