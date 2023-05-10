#include "barrier_walk.hpp"

class DikinWalk: public BarrierWalk{

    public:

        DikinWalk() : BarrierWalk(){}

        void initialize(MatrixXd A_p, VectorXd b_p, float r){
            float constant = (r * r)/b_p.rows();
            float td = (-0.5 / constant);
            float ts = sqrt(constant);
            A = A_p;
            b = b_p;
            setTs(ts);
            setTd(td);
        }
        void generateWeight(VectorXd& x);

        void printType();
};