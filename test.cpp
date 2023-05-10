#include "dikin_walk.hpp"
#include "dikinls_walk.hpp"
#include "vaidya_walk.hpp"
#include "john_walk.hpp"
#include "ball_walk.hpp"
#include "hitrun_walk.hpp"

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
}