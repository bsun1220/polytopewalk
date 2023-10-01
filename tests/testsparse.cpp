#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cstring>
#include <sparserandomwalk/PolytopeWalk.hpp>

using namespace std;

static string one() {
    return "1";
}

TEST_CASE( "Assert that something is true (pass)", "[require]" ) {
    REQUIRE( one() == "1" );
}

TEST_CASE( "Check Ball Walker Results (pass)", "[require]" ) {
    MatrixXd A (4,2);
    A << 1, 0, -1, 0, 0, 1, 0, -1;
    SparseMatrixXd AS (4,2);
    AS = A.sparseView();

    VectorXd b(4);
    b << 1, 0, 1, 0;

    VectorXd init(2);
    init << 0.5, 0.5;
    const int steps = 1000;

    BallWalk ball(0.4);
    MatrixXd res_ball = ball.generateCompleteWalk(steps, init, AS, b);

    REQUIRE(res_ball.rows() == steps);
    REQUIRE(res_ball.cols() == 2);
    REQUIRE(res_ball.minCoeff() > 0);
    REQUIRE(res_ball.maxCoeff() < 1);

}

TEST_CASE( "Check Dikin Walker Results (pass)", "[require]" ) {
    MatrixXd A (4,2);
    A << 1, 0, -1, 0, 0, 1, 0, -1;
    SparseMatrixXd AS (4,2);
    AS = A.sparseView();

    VectorXd b(4);
    b << 1, 0, 1, 0;

    VectorXd init(2);
    init << 0.5, 0.5;
    const int steps = 1000;

    DikinWalk dikin (0.5);
    MatrixXd res_dikin = dikin.generateCompleteWalk(steps, init, AS, b);

    REQUIRE(res_dikin.rows() == steps);
    REQUIRE(res_dikin.cols() == 2);
    REQUIRE(res_dikin.minCoeff() > 0);
    REQUIRE(res_dikin.maxCoeff() < 1);
}

TEST_CASE( "Check Vaidya Walker Results (pass)", "[require]" ) {
    MatrixXd A (4,2);
    A << 1, 0, -1, 0, 0, 1, 0, -1;
    SparseMatrixXd AS (4,2);
    AS = A.sparseView();

    VectorXd b(4);
    b << 1, 0, 1, 0;

    VectorXd init(2);
    init << 0.5, 0.5;
    const int steps = 1000;

    VaidyaWalk vaidya(0.5);
    MatrixXd res_vaidya = vaidya.generateCompleteWalk(steps, init, AS, b);

    REQUIRE(res_vaidya.rows() == steps);
    REQUIRE(res_vaidya.cols() == 2);
    REQUIRE(res_vaidya.minCoeff() > 0);
    REQUIRE(res_vaidya.maxCoeff() < 1);

}


TEST_CASE( "Check HitRun Walker Results (pass)", "[require]" ) {
    MatrixXd A (4,2);
    A << 1, 0, -1, 0, 0, 1, 0, -1;
    SparseMatrixXd AS (4,2);
    AS = A.sparseView();

    VectorXd b(4);
    b << 1, 0, 1, 0;

    VectorXd init(2);
    init << 0.5, 0.5;
    const int steps = 1000;

    HitAndRunWalk hr(0.01, 0.1);
    MatrixXd res_hr = hr.generateCompleteWalk(steps, init, AS, b);

    REQUIRE(res_hr.rows() == steps);
    REQUIRE(res_hr.cols() == 2);
    REQUIRE(res_hr.minCoeff() > 0);
    REQUIRE(res_hr.maxCoeff() < 1);
}

TEST_CASE( "Check DikinLS Walker Results (pass)", "[require]" ) {
    MatrixXd A (4,2);
    A << 1, 0, -1, 0, 0, 1, 0, -1;
    SparseMatrixXd AS = A.sparseView();

    VectorXd b(4);
    b << 1, 0, 1, 0;

    VectorXd init(2);
    init << 0.5, 0.5;
    const int steps = 1000;

    DikinLSWalk dikinls (0.1, 100, 0.001, 0.4);
    MatrixXd res_dikinls = dikinls.generateCompleteWalk(steps, init, AS, b);

    REQUIRE(res_dikinls.rows() == steps);
    REQUIRE(res_dikinls.cols() == 2);
    REQUIRE(res_dikinls.minCoeff() > 0);
    REQUIRE(res_dikinls.maxCoeff() < 1);

}


TEST_CASE( "Check John Walker Results (pass)", "[require]" ) {
    MatrixXd A (4,2);
    A << 1, 0, -1, 0, 0, 1, 0, -1;
    SparseMatrixXd AS = A.sparseView();

    VectorXd b(4);
    b << 1, 0, 1, 0;

    VectorXd init(2);
    init << 0.5, 0.5;
    const int steps = 1000;

    JohnWalk john(0.1, 100, 0.001, 0.4);
    MatrixXd res_john = john.generateCompleteWalk(steps, init, AS, b);

    REQUIRE(res_john.rows() == steps);
    REQUIRE(res_john.cols() == 2);
    REQUIRE(res_john.minCoeff() > 0);
    REQUIRE(res_john.maxCoeff() < 1);

}

TEST_CASE( "Check CentralPointFinder Results (pass)", "[require]" ) {
    CentralPointFinder cpf (10000, 0.0001, 10000, 0.0001, 0.0001);
    
    MatrixXd A1 (4,2);
    A1 << 1, 0, -1, 0, 0, 1, 0, -1;

    VectorXd b1(4);
    b1 << 1, 0, 1, 0;

    VectorXd x1 = cpf.getInitialPoint(A1, b1);

    REQUIRE(x1.rows() == 2);
    REQUIRE(x1(0) == Catch::Approx(0.5).epsilon(0.01));
    REQUIRE(x1(1) == Catch::Approx(0.5).epsilon(0.01));

}

TEST_CASE( "Check FacialReduction Results (pass)", "[require]" ) {
    FacialReduction fr; 

    MatrixXd A1 (6, 3);
    A1 << 1, 1, 0, -1, -1, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;

    VectorXd b1(6);
    b1 << 1, -1, 1, 1, 1, 1;

    problem_result res1 = fr.reduce(A1, b1);

    REQUIRE(res1.reduced == true);
    REQUIRE(res1.reduced_A.cols() == 2);
    REQUIRE(res1.reduced_A.rows() == 4);
    REQUIRE(res1.reduced_b.rows() == 4);

    MatrixXd A2(6,3);
    A2 << 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;

    VectorXd b2(6);
    b2 << 1, 1, 0, 0, 0, 0;

    problem_result res2 = fr.reduce(A2, b2);
    
    REQUIRE(res2.reduced == true);
    REQUIRE(res2.reduced_A.cols() == 1);
    REQUIRE(res2.reduced_A.rows() == 2);
    REQUIRE(res2.reduced_b.rows() == 2);

    
    MatrixXd A3(4,2);
    A3 << 1, 0, -1, 0, 0, 1, 0, -1;

    VectorXd b3(4);
    b3 << 1, 0, 1, 0;
    problem_result res3 = fr.reduce(A3, b3);

    REQUIRE(res3.reduced == false);

}


TEST_CASE( "Check PolytopeWalk Results (pass)", "[require]" ){
    MatrixXd A (6, 3);
    A << 1, 1, 0, -1, -1, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;
    SparseMatrixXd AS (6,3);
    AS = A.sparseView();

    VectorXd b(6);
    b << 1, -1, 1, 1, 1, 1;

    DikinWalk dikin (0.5);
    FacialReduction fr;
    CentralPointFinder cpf (10000, 0.0001, 10000, 0.0001, 0.0001);
    MatrixXd res = fullWalkRun(AS, b, 1000, &dikin, &fr, &cpf);

    REQUIRE(res.rows() == 1000);
    REQUIRE(res.cols() == 3);

}