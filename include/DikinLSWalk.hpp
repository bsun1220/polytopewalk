#include "BarrierWalk.hpp"

class DikinLSWalk: public BarrierWalk{

    public:
        /**
         * @brief constructor for weighted Dikin Walk class
         * @param ss 
         * @param mi 
         * @param gl 
         * @param rp
         */
        DikinLSWalk(const float ss, const int mi, const float gl, const float rp) : STEPSIZE(ss), MAXITER(mi), GRADLIM(gl), BarrierWalk(rp){
            
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
        const float STEPSIZE;

        /**
         * @brief max number of iterations in gradient descent
         */
        const int MAXITER;

        /**
         * @brief stops gradient descent if it reaches under this number
         */
        const float GRADLIM;

        /**
         * @brief generate weights when calculating Hessian matrix
         * @param x point in polytope to generate DikinLS weight
         * @param A polytope matrix
         * @param b polytope vector
         * @return void (update global variable weights)
         */
        void generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b);
};