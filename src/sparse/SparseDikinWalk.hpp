#ifndef CONSDIKIN_HPP
#define CONSDIKIN_HPP

#include "SparseBarrierWalk.hpp"

class SparseDikinWalk : public SparseBarrierWalk{

    public:
        /**
         * @brief constructor for Dikin Walk class
         * @param err error constant
         * @param r spread parameter
         * @param thin thin parameter
         */
        SparseDikinWalk(double r, double err = 1e-6, int thin = 1) : SparseBarrierWalk(r, err, thin) {}

        /**
         * @brief generate weight (identity matrix)
         * @param x slack variable
         * @param A polytope constraint
         * @param k k values >= 0 constraint
         * @return SparseMatrixXd
         */
        SparseMatrixXd generateWeight(
            const VectorXd& x, 
            const SparseMatrixXd& A,
            int k
        ) override; 
    
    protected:

        /**
         * @brief Distribution constant
         * @param d polytope matrix 
         * @param n polytope vector
         * @return void
         */
        void setDistTerm(int d, int n) override;

};

#endif

