#include "BarrierWalk.hpp"

class DikinLSWalk: public BarrierWalk{

    public:
        DikinLSWalk(float ss, float mi, float gl) : BarrierWalk(){
            step_size = ss;
            max_iter = mi;
            grad_lim = gl;
        }

        void initialize(MatrixXd A_p, VectorXd b_p, float r){
            float constant = (r * r)/b_p.rows();
            float td = (-0.5 / constant);
            float ts = sqrt(constant);
            A = A_p;
            b = b_p;
            BarrierWalk::setTs(ts);
            BarrierWalk::setTd(td);
        }
        void generateWeight(VectorXd& x);
        void printType();
    
    protected:
        float step_size {};
        float max_iter {};
        float grad_lim {};
        void gradientDescent(VectorXd& x, float adj, int sim, float grad_lim);
};