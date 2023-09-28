
#ifndef DIKINWALK_HPP
#define DIKINWALK_HPP

#include "BarrierWalk.hpp"

class DikinWalk: public BarrierWalk{

    public:

        DikinWalk(const float rp) : BarrierWalk(rp){}

        /**
         * @brief returns weight for DikinWalk (Identity Matrix)
         * @param x point
         * @param A polytope matrix
         * @param b polytope matrix
         * @return void
         */
        void generateWeight(const VectorXd& x, const SparseMatrixXd& A, const VectorXd&b) override;

        /**
         * @brief print dikin
         * @return void
         */
        void printType() override;
};


#endif