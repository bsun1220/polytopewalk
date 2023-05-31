#include "Reducer.hpp"

/**
 * @brief object for whether a solution was found or not and its corresponding value
 */
struct z_result{
    bool found_sol; 
    VectorXd z; 
};

/**
 * @brief object for edited versions of A and b during facial reduction process
 */
struct fr_result{
    MatrixXd A;
    VectorXd b; 
};

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
        MatrixXd equalConversion(MatrixXd A);

         /**
         * @brief finds a vector z satisfying A^Ty = [0 z], z in R^n, z >= 0, z != 0, <b, y> = 0
         * @param A
         * @return Matrix
         */
        z_result findZ(MatrixXd newA, VectorXd b, int x_dim);

         /**
         * @brief Finds a Matrix V to convert AVv = b after receiving v
         * @param z
         * @return Matrix
         */
        MatrixXd facialReduction(VectorXd z);

         /**
         * @brief iteratively reduces dimension of the problem using recursion
         * @param A
         * @param b
         * @return Matrix A, Vector b
         */
        fr_result entireFacialReductionStep(MatrixXd A, VectorXd b, int x_dim);

         /**
         * @brief completes the square to find new A, b for the reduced polytope
         * @param A
         * @return Matrix A, Vector b
         */
        fr_result reduceSampling(MatrixXd M, VectorXd b, int delta_dim);

};