#include "BarrierWalk.hpp"

class JohnWalk: public BarrierWalk{

    public:
        /**
         * @brief constructor for weighted John walk class
         * @param ss 
         * @param mi 
         * @param gl 
         */
        JohnWalk(float ss, float mi, float gl) : BarrierWalk(){
            step_size = ss;
            max_iter = mi;
            grad_lim = gl;
        }

        /**
         * @brief Initialize values (because prior to Reducer, it is unknown what these values are)
         * @param A_p (Matrix for polytope) {x | Ax <= b}
         * @param b_p (Vector for polytope) {x | Ax <= b}
         * @param r_p general indicator of spread
         * @return void
         */
        void initialize(MatrixXd A_p, VectorXd b_p, float r){
            float constant = (r * r)/b_p.rows();
            float td = (-0.5 / constant);
            float ts = sqrt(constant);
            A = A_p;
            b = b_p;
            BarrierWalk::setTs(ts);
            BarrierWalk::setTd(td);
        }

        /**
         * @brief print dikinls
         * @return void
         */
        void printType();
    
    protected:
        /**
         * @brief step size for gradient descent
         */
        float step_size {};

        /**
         * @brief max number of iterations in gradient descent
         */
        float max_iter {};

        /**
         * @brief stops gradient descent if it reaches under this number
         */
        float grad_lim {};

        /**
         * @brief generate weights when calculating Hessian matrix
         * @param x point in polytope to generate DikinLS weight
         * @return void (update global variable weights)
         */
        void generateWeight(VectorXd& x);

         /**
         * @brief solves an optimization problem to generate weight x
         */
        void gradientDescent(VectorXd& x, float adj, int sim, float grad_lim);
};