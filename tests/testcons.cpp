#include "constraintwalk/SparseFacialReduction.hpp"
#include "constraintwalk/SparseCenter.hpp"
#include "constraintwalk/ConstraintDikinWalk.hpp"
#include "constraintwalk/ConstraintVaidyaWalk.hpp"
#include "constraintwalk/ConstraintDikinLSWalk.hpp"
#include "constraintwalk/ConstraintJohnWalk.hpp"

int main(){

    SparseMatrixXd A (4, 6);
    A.coeffRef(0, 0) = 1;
    A.coeffRef(0, 2) = 1;
    A.coeffRef(1, 1) = 1;
    A.coeffRef(1, 3) = 1;
    A.coeffRef(2, 0) = -1;
    A.coeffRef(2, 4) = 1;
    A.coeffRef(3, 1) = -1;
    A.coeffRef(3, 5) = 1;
    VectorXd b (4);
    b << 1, 1, 1, 1;
    VectorXd x (6);
    x << 0, 0, 1, 1, 1, 1;

    ConstraintDikinLSWalk cd(0.0001, 0.4, 0.01, 0.1, 1000);
    MatrixXd res = cd.generateCompleteWalk(10, x, A, b, 4);
    cout << res << endl;
}
