#ifndef INIT_H
#define INIT_H
#include "Common.hpp"

class Initializer{
    public:

    Initializer (){};

    virtual VectorXd getInitialPoint(MatrixXd A, VectorXd b);
};


#endif


