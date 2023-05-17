#include "BarrierWalk.hpp"

class DikinWalk: public BarrierWalk{

    public:

        DikinWalk(float rp) : BarrierWalk(rp){}

        /**
         * @brief returns weight for DikinWalk (Identity Matrix)
         * @param x point
         * @param A polytope matrix
         * @param b polytope matrix
         * @return void
         */
        void generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd&b);

        /**
         * @brief print dikin
         * @return void
         */
        void printType();
};