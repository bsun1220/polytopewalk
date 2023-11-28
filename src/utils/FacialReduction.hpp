
#ifndef FACIALREDUCTION_HPP
#define FACIALREDUCTION_HPP

#include "Reducer.hpp"

class FacialReduction: public Reducer{

    public:
        FacialReduction() : Reducer(){};

        /**
         * @brief converts A for the problem
         * @param A
         * @param b
         * @return problem_result object for reducedA, reducedb, reduced value, originalA, originalb
         */
        problem_result reduce(MatrixXd A, VectorXd b) override;

    protected:
        /**
         * @brief converts A for the problem
         * @param A
         * @return Matrix
         */
        MatrixXd equalConversion(const MatrixXd& A);

         /**
         * @brief finds a vector z satisfying A^Ty = [0 z], z in R^n, z >= 0, z != 0, <b, y> = 0
         * @param A
         * @return Matrix
         */
        z_result findZ(const MatrixXd& newA, const VectorXd& b, int x_dim);

         /**
         * @brief Finds a Matrix V to convert AVv = b after receiving v
         * @param z
         * @param x_dim
         * @return Matrix
         */
        MatrixXd pickV(const VectorXd& z, int x_dim);

         /**
         * @brief Finds the Projection Matrix
         * @param AV
         * @return Matrix
         */
        MatrixXd pickP(const MatrixXd& AV);

         /**
         * @brief iteratively reduces dimension of the problem using recursion
         * @param A
         * @param b
         * @return Matrix A, Vector b
         */
        fr_result entireFacialReductionStep(MatrixXd A, VectorXd b, int x_dim);

};

#endif