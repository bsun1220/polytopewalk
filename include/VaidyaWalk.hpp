#include "BarrierWalk.hpp"

class VaidyaWalk: public BarrierWalk{

    public:

        VaidyaWalk() : BarrierWalk(){}  

        /**
         * @brief Initialize values (because prior to Reducer, it is unknown what these values are)
         * @param A_p (Matrix for polytope) Ax <= b
         * @param b_p (Vector for polytope) Ax <= b
         * @param r_p general indicator of spread
         * @return void
         */
        void initialize(MatrixXd A_p, VectorXd b_p, float r){
            float constant = (r * r)/sqrt(A_p.cols() * A_p.rows());
            float td = (-0.5 / constant);
            float ts = sqrt(constant);
            A = A_p;
            b = b_p;
            BarrierWalk::setTs(ts);
            BarrierWalk::setTd(td);
        }

        /**
         * @brief print general type 
         * @return void
         */
        void printType();
    
    protected:
         /**
         * @brief global variable to update dikin hessian
         */
        MatrixXd dhess {};

         /**
         * @brief returns weight for Vaidya Walk
         * @param x
         * @return void
         */
        void generateWeight(VectorXd& x);

        /**
         * @brief returns unweight for Dikin Hessian around vector x
         * @param x
         * @return void
         */
        void generateDikinHessian(VectorXd& x);
};