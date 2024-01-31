#ifndef CONSDIKINLS_HPP
#define CONSDIKINLS_HPP

#include "ConstraintBarrierWalk.hpp"
#include "LeverageScore.hpp"

class ConstraintDikinLSWalk : public ConstraintBarrierWalk{

    public:
        ConstraintDikinLSWalk(double err, double r, double g_lim, double step_size, int max_iter) : G_LIM(g_lim), STEP_SIZE(step_size), MAX_ITER(max_iter), ConstraintBarrierWalk(err, r) {}
    
    protected:
        SparseMatrixXd generateG(
            const VectorXd& x, 
            const SparseMatrixXd& A,
            int k
        ) override; 

        void setDistTerm(int d, int n) override;
        const double G_LIM;
        const double STEP_SIZE;
        const int MAX_ITER;

};

#endif