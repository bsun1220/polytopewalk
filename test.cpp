#include "PolytopeWalk.hpp"

int main(){

    MatrixXd A (4,2);
    A << 1, 0, -1, 0, 1, 1, -1, -1;

    VectorXd b(4);
    b << 1, 1, 1, -1;

    //BallWalk walk;
    JohnWalk walk(0.01, 100, 0.01);
    FacialReduction fr;

    CentralPointFinder cpf(10000, 0.0001, 10000, 0.0001, 0.001);

    MatrixXd res = fullWalkRun(A, b, 0.4, 100, walk, fr, cpf);
    cout << res << endl;

    /*
    MatrixXd A (4,2);
    A << 1, 0, -1, 0, 0, 1, 0, -1;
    VectorXd b (4);
    b << 1, 0, 1, 0;
    VectorXd x (2);
    x << 0.5, 0.5;


    MatrixXd A (6,3);
    A << 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;
    VectorXd b(6);
    b << 1, 0, 1, 0, 1, 0;
    VectorXd x (3);
    x << 0.5, 0.5, 0.5;

    DikinWalk d;
    d.initialize(A, b, 0.4);
    cout <<  d.generateCompleteWalk(10, x) << endl;

    DikinLSWalk dl (0.1, 100, 0.01);
    dl.initialize(A, b, 0.4);
    cout <<  dl.generateCompleteWalk(10, x) << endl;

    VaidyaWalk v;
    v.initialize(A, b, 0.4);
    cout << v.generateCompleteWalk(10, x) << endl;

    JohnWalk j (0.1, 100, 0.01);
    j.initialize(A, b, 0.4);
    cout <<  j.generateCompleteWalk(10, x) << endl;

    BallWalk ba;
    ba.initialize(A, b, 0.1);
    cout << ba.generateCompleteWalk(10, x) << endl;

    HitAndRunWalk hr (0.2);
    hr.initialize(A, b, 0.1);
    cout << hr.generateCompleteWalk(10, x) << endl;

    CentralPointFinder cpf (0.0001, 0.0001, 10000.0, 0.00001, 0.01);
    cout << cpf.getInitialPoint(A, b) << endl;

    MatrixXd newA (6, 3);
    newA << 1, 1, 0, -1, -1, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;

    VectorXd newB (6);
    newB << 1, -1, 1, 1, 1, 1;

    FacialReduction fr;
    problem_result res = fr.reduce(newA, newB);

    cout << res.reduced_A << endl;
    cout << res.reduced_b << endl;

    MatrixXd newA2 (6,3);
    newA2 << 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;

    VectorXd newB2(6);
    newB2 << 1, 1, 0, 0, 0, 0;

    res = fr.reduce(newA2, newB2);
    cout << res.reduced_A << endl;
    cout << res.reduced_b << endl;

    MatrixXd newA3 (4,2);
    newA3 << 1, 0, -1, 0, 1, 1, -1, -1;

    VectorXd newB3(4);
    newB3 << 1, 1, 1, -1;

    res = fr.reduce(newA3, newB3);
    cout << res.reduced_A << endl;
    cout << res.reduced_b << endl;
    */

    
}


