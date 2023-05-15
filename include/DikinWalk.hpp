#include "BarrierWalk.hpp"

class DikinWalk: public BarrierWalk{

    public:

        DikinWalk() : BarrierWalk(){}

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
            setTs(ts);
            setTd(td);
        }
        /**
         * @brief returns weight for DikinWalk (Identity Matrix)
         * @return void
         */
        void generateWeight(VectorXd& x);

        /**
         * @brief print dikin
         * @return void
         */
        void printType();
};