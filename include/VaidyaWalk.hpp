#include "BarrierWalk.hpp"

class VaidyaWalk: public BarrierWalk{

    public:

        VaidyaWalk(const float rp) : BarrierWalk(rp){}  

        /**
         * @brief print general type 
         * @return void
         */
        void printType() override;
    
    protected:
         /**
         * @brief global variable to update dikin hessian
         */
        MatrixXd dhess {};

         /**
         * @brief returns weight for Vaidya Walk
         * @param x
         * @param A polytope matrix
         * @param b polytope vector
         * @return void
         */
        void generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b) override;

        /**
         * @brief returns unweight for Dikin Hessian around vector x
         * @param x
         * @param A
         * @param b
         * @return void
         */
        void generateDikinHessian(const VectorXd& x, const MatrixXd& A, const VectorXd& b);
};