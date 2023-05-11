#include "DikinWalk.hpp"
#include "DikinLSWalk.hpp"
#include "VaidyaWalk.hpp"
#include "JohnWalk.hpp"
#include "BallWalk.hpp"
#include "HitRunWalk.hpp"

#include "CentralPointFinder.hpp"

int main(){
    /*
    MatrixXd A (4,2);
    A << 1, 0, -1, 0, 0, 1, 0, -1;
    VectorXd b (4);
    b << 1, 0, 1, 0;
    VectorXd x (2);
    x << 0.5, 0.5;*/

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
}