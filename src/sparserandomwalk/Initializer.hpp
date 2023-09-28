#ifndef INITIALIZER_HPP
#define INITIALIZER_HPP
#include "Common.hpp"

class Initializer{
    public:

    Initializer (){};
    /**
     * @brief template get the initial point to start polytope walk
     * @param A
     * @param b
     * @return Vector
     */
    virtual VectorXd getInitialPoint(MatrixXd A, VectorXd b);
};


#endif


