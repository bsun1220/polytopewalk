
#include "utils/FacialReduction.hpp"

int main(){

    FacialReduction fr;

    MatrixXd A1 (6, 3);
    A1 << 1, 1, 0, -1, -1, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;

    VectorXd b1(6);
    b1 << 1, -1, 1, 1, 1, 1;

    problem_result res1 = fr.reduce(A1, b1);

    MatrixXd A2(6,3);
    A2 << 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;

    VectorXd b2(6);
    b2 << 1, 1, 0, 0, 0, 0;

    problem_result res2 = fr.reduce(A2, b2);

    MatrixXd A3(4,2);
    A3 << 1, 0, -1, 0, 0, 1, 0, -1;

    VectorXd b3(4);
    b3 << 1, 0, 1, 0;
    problem_result res3 = fr.reduce(A3, b3);
}
