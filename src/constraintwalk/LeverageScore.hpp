#ifndef LEVSCORE_HPP
#define LEVSCORE_HPP
#include "Common.hpp"

class LeverageScore{
    public:

    LeverageScore (){};
    /**
     * @brief Get the Leverage Score approximate calculation
     * @param A
     * @param W
     * @param b
     * @return Vector
     */
    VectorXd generate(const SparseMatrixXd& A, const SparseMatrixXd& W, const VectorXd& x);
};

#endif