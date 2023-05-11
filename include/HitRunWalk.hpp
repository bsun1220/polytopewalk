#include "RandomWalk.hpp"

class HitAndRunWalk: public RandomWalk{

    float err {};
    public:
        HitAndRunWalk(float err_p) : RandomWalk() {
            err = err_p; 
        }

        void initialize(MatrixXd A_p, VectorXd b_p, float r_p){
            A = A_p;
            b = b_p;
            r = r_p;
        }

        double distance(VectorXd& x, VectorXd&y);

        double binarySearch(VectorXd direction, VectorXd x);

        MatrixXd generateCompleteWalk(int num_steps, VectorXd x);

};