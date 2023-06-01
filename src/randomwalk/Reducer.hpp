#ifndef REDUCER_HPP
#define REDUCER_HPP
#include "Common.hpp"
#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>
using namespace ifopt;

struct problem_result{
    MatrixXd reduced_A;
    VectorXd reduced_b; 
    bool reduced;
    VectorXd b_tilde;
    MatrixXd M;
};

class Reducer{
    public:

    Reducer (){};

    MatrixXd makeFullRank(MatrixXd mat);
    virtual problem_result reduce(MatrixXd A, VectorXd b);

};


#endif
