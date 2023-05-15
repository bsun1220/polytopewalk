#include "RandomWalk.hpp"

class HitAndRunWalk: public RandomWalk{

    public:
        /**
         * @brief Hit and Run implementation constructor
         * 
         */
        HitAndRunWalk(float err_p) : RandomWalk() {
            err = err_p; 
        }

        /**
         * @brief Initialize values (because prior to Reducer, it is unknown what these values are)
         * @param A_p (Matrix for polytope) {x | Ax <= b}
         * @param b_p (Vector for polytope) {x | Ax <= b}
         * @param r_p general indicator of spread
         * @return void
         */
        void initialize(MatrixXd A_p, VectorXd b_p, float r_p){
            A = A_p;
            b = b_p;
            r = r_p;
        }
        /**
         * @brief Generate values from the walk
         * @param num_steps number of steps wanted to take
         * @param x initial starting point
         * @return Matrix
         */
        MatrixXd generateCompleteWalk(const int num_steps, VectorXd& x);

         /**
         * @brief print general type 
         * @return void
         */
        void printType();
    
    protected:
        /**
         * @brief relative error of the binary search operation
         */
        float err {};

        /**
         * @brief get distance between vectors x and y
         * @param x
         * @param y
         * @return double 
         */
        double distance(VectorXd& x, VectorXd&y);

        /**
         * @brief runs binary search to find a suitable chord intersection with the polytope
         * @param direction (random direction variable)
         * @param x (starting point)
         * @return double 
         */
        double binarySearch(VectorXd direction, VectorXd& x);

};