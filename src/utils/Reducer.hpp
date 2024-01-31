#ifndef REDUCER_HPP
#define REDUCER_HPP
#include "Common.hpp"
#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>
using namespace ifopt;

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

struct problem_result{
    MatrixXd reduced_A;
    VectorXd reduced_b; 
    bool reduced;
    VectorXd z1;
    MatrixXd Q;
};

class Reducer{
    public:

    Reducer (){};

    MatrixXd makeFullRank(const MatrixXd& mat);
    virtual problem_result reduce(MatrixXd A, VectorXd b);

};
class ExVariables: public VariableSet{
    public:
    VectorXd x;

    ExVariables(int num_dim, string name, const VectorXd& init) : VariableSet(num_dim, name){
        x = init;
    } 
    void SetVariables(const VectorXd& val) override
    {
        x = val;
    }

    VectorXd GetValues() const override{
        return x;
    }

    VecBound GetBounds() const override{
        VecBound bounds(GetRows());
        for(int i = 0; i < x.rows(); i++){
            bounds.at(i) = NoBound;
        }
        return bounds;
    }

};

class ExConstraint1 : public ConstraintSet{
    public:
    MatrixXd A;
    ExConstraint1(int num_dim, string name, const MatrixXd& A_param) : ConstraintSet(num_dim, name){
        // A is d by n matrix
        A = A_param;
    }

    VectorXd GetValues() const override{
        VectorXd input = GetVariables()->GetComponent("var_set1")->GetValues();
        return A * input; 
    }

    VecBound GetBounds() const override{
        VecBound b(GetRows());
        for(int i = 0; i < A.rows(); i++){
            b.at(i) = Bounds(-inf, 0);
        }
        return b;
    }
    void FillJacobianBlock (string var_set, Jacobian& jac_block) const override{

        for(int i = 0; i < A.rows(); i++){
            for(int j = 0; j < A.cols(); j++){
                if (A.coeff(i, j) != 0){
                    if (jac_block.coeffRef(i, j) != 0){
                        return;
                    }
                    jac_block.coeffRef(i, j) = A.coeff(i, j);
                }
            }
        }
    }

};

class ExConstraint2 : public ConstraintSet{
    public:
    MatrixXd A;
    VectorXd b;
    ExConstraint2(int num_dim, string name, const MatrixXd& A_param, const VectorXd& b_param) : ConstraintSet(num_dim, name){
        A = A_param;
        b = b_param;

    }

    VectorXd GetValues() const override{
        VectorXd input = GetVariables()->GetComponent("var_set1")->GetValues();
        return A * input;
    }

    VecBound GetBounds() const override{
        VecBound bound(GetRows());
        for(int i = 0; i < b.rows(); i++){
            bound.at(i) = Bounds(b(i),b(i));
        }
        return bound;
    }

    void FillJacobianBlock (string var_set, Jacobian& jac_block) const override{
        for(int i = 0; i < A.rows(); i++){
            for(int j = 0; j < A.cols(); j++){
                if (A.coeff(i, j) != 0){
                    if (jac_block.coeffRef(i, j) != 0){
                        return;
                    }
                    jac_block.coeffRef(i, j) = A.coeff(i, j);
                }
            }
        }
    }

};

class ExCost : public CostTerm{
    public:

    ExCost(string name) : CostTerm(name) {

    }
    
    double GetCost() const override
    {
        VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();
        return x(x.rows() - 1);
    };
    void FillJacobianBlock (string var_set, Jacobian& jac_block) const override{
        VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();
        jac_block.coeffRef(0, x.rows() - 1) = 1; 
    }
};


#endif
