#include "Reducer.hpp"

struct z_result{
    bool found_sol; 
    VectorXd z; 
};

struct fr_result{
    MatrixXd A;
    VectorXd b; 
};

class FacialReduction: public Reducer{

    public:
        FacialReduction() : Reducer(){};
        problem_result reduce(MatrixXd A, VectorXd b);


    protected:
        MatrixXd equalConversion(MatrixXd A);
        z_result findZ(MatrixXd newA, VectorXd b, int x_dim);
        MatrixXd facialReduction(VectorXd z);
        fr_result entireFacialReductionStep(MatrixXd A, VectorXd b, int x_dim);
        fr_result reduceSampling(MatrixXd M, VectorXd b, int delta_dim);


};