
#include "Initializer.hpp"


class CentralPointFinder : public Initializer {
    public:

    /**
    * @brief Initializer for an algorithm which finds central point of a polytope
    * @param ss initial starting point for lp solver
    * @param ts initial step size
    * @param te final step size
    * @param err_term fixed term to prevent divide by 0 error
    * @param grad_lim run gradient descent until we reach this number
    */
    CentralPointFinder(const double ss, const double ts, const double te, 
                       const double err_term, const double grad_lim) : 
                       ss_(ss), ts_(ts), te_(te), err_term_(err_term), grad_lim_(grad_lim), Initializer() {}


    /**
    * @brief let tilde(a_i) = [a_i .. -1], y = [x ... delta], and c = [0... 0, 1]
    * we solve the optimization problem min c^Ty s.t. tilde(a_i)^Ty - b <= 0
    * @param A Matrix in polytope {x | Ax <= b}
    * @param b Vector in polytope {x | Ax <= b}
    * @return Vector (The geometrical center of the polytope)
    */
    VectorXd getInitialPoint(MatrixXd A, VectorXd b);


    protected:
        /**
        * @brief starting error point in lp
        **/
        const double ss_;

        /**
        * @brief starting initial iterator
        **/
        const double ts_;

        /**
        * @brief final iterator 
        **/
        const double te_;
        /**
        * @brief error term to prevent divide 0 error
        **/
        const double err_term_;

        /**
        * @brief running gradient descent until this limit
        **/
        const double grad_lim_;

};