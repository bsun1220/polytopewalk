#ifndef CPF_HPP
#define CPF_HPP

#include "Common.hpp"


class DenseCenter {
    public:

    DenseCenter(int max_iter = 3000, double tol = 1e-8, double s_max = 100) : MAX_ITER(max_iter), TOL(tol), S_MAX(s_max) {};

    /**
     * @brief Finds analytical center Ax <= b
     * @param A polytope matrix (Ax <= b)
     * @param b polytope vector (Ax <= b)
     * @return VectorXd 
     */
    VectorXd getInitialPoint(MatrixXd A, VectorXd b);

    protected:
        const int MAX_ITER;
        const double TOL;
        const double S_MAX;

};

#endif