#include "BarrierWalk.hpp"

class VaidyaWalk: public BarrierWalk{

    public:

        VaidyaWalk() : BarrierWalk(){}

        void initialize(MatrixXd A_p, VectorXd b_p, float r){
            float constant = (r * r)/sqrt(A_p.cols() * A_p.rows());
            float td = (-0.5 / constant);
            float ts = sqrt(constant);
            A = A_p;
            b = b_p;
            BarrierWalk::setTs(ts);
            BarrierWalk::setTd(td);
        }
        void generateWeight(VectorXd& x);
        void generateDikinHessian(VectorXd& x);
        void printType();
    
    protected:
        MatrixXd dhess {};
};