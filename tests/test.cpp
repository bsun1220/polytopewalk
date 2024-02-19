#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "utils/FullWalkRun.hpp"
#include <cstring>

struct sparse_polytope{
    SparseMatrixXd A;
    VectorXd b; 
    int k;
};

sparse_polytope generate_simplex(){
    SparseMatrixXd simplex_A (1, 3);
    simplex_A.coeffRef(0, 0) = 1;
    simplex_A.coeffRef(0, 1) = 1;
    simplex_A.coeffRef(0, 2) = 1;
    VectorXd simplex_b (1);
    simplex_b << 1;
    sparse_polytope result; 
    result.A = simplex_A;
    result.b = simplex_b;
    result.k = 3; 
    return result;
}

sparse_polytope generate_hc(){
    SparseMatrixXd hc_A (4, 6);
    hc_A.coeffRef(0, 0) = 1;
    hc_A.coeffRef(0, 2) = 1;
    hc_A.coeffRef(1, 1) = 1;
    hc_A.coeffRef(1, 3) = 1;
    hc_A.coeffRef(2, 0) = -1;
    hc_A.coeffRef(2, 4) = 1;
    hc_A.coeffRef(3, 1) = -1;
    hc_A.coeffRef(3, 5) = 1; 

    VectorXd hc_b (4);
    hc_b << 1, 1, 1, 1;
    sparse_polytope result; 
    result.A = hc_A;
    result.b = hc_b;
    result.k = 4; 
    return result;
}

sparse_polytope generate_birkhoff(){
    SparseMatrixXd birk_A (3, 4);
    birk_A.coeffRef(0, 0) = 1;
    birk_A.coeffRef(0, 1) = 1;
    birk_A.coeffRef(1, 2) = 1;
    birk_A.coeffRef(1, 3) = 1;
    birk_A.coeffRef(2, 0) = 1;
    birk_A.coeffRef(2, 2) = 1; 

    VectorXd birk_b (3);
    birk_b << 1, 1, 1;
    sparse_polytope result; 
    result.A = birk_A;
    result.b = birk_b;
    result.k = 4; 
    return result;
}

sparse_polytope simplex = generate_simplex();
sparse_polytope hc = generate_hc();
sparse_polytope birk = generate_birkhoff();

TEST_CASE( "Check Facial Reduction Algorithm", "[require]" ) {
    FacialReduction fr;
    res simplex_dense = fr.reduce(simplex.A, simplex.b, simplex.k, false);
    res hc_dense = fr.reduce(hc.A, hc.b, hc.k, false);
    res birk_dense = fr.reduce(birk.A, birk.b, birk.k, false);

    REQUIRE(((simplex_dense.dense_A.rows() == 3) && (simplex_dense.dense_A.cols() == 2)));
    REQUIRE(simplex_dense.dense_b.rows() == 3);
    REQUIRE(((hc_dense.dense_A.rows() == 4) && (hc_dense.dense_A.cols() == 2)));
    REQUIRE(hc_dense.dense_b.rows() == 4);
    REQUIRE(((birk_dense.dense_A.rows() == 4) && (birk_dense.dense_A.cols() == 1)));
    REQUIRE(birk_dense.dense_b.rows() == 4);

    MatrixXd A1 (6, 3);
    A1 << 1, 1, 0, -1, -1, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;
    MatrixXd temp (6, 9);
    temp << A1, VectorXd::Ones(6).asDiagonal().toDenseMatrix();
    SparseMatrixXd SA1 = temp.sparseView();
    VectorXd b1(6);
    b1 << 1, -1, 1, 1, 1, 1;

    res test1a = fr.reduce(SA1, b1, 6, true);
    REQUIRE((test1a.sparse_A.rows() == 5 && test1a.sparse_A.cols() == 7));
    res test1b = fr.reduce(SA1, b1, 6, false);
    REQUIRE((test1b.dense_A.rows() == 4 && test1b.dense_A.cols() == 2));

    MatrixXd A2(6,3);
    A2 << 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;
    MatrixXd temp2 (6, 9);
    temp2 << A2, VectorXd::Ones(6).asDiagonal().toDenseMatrix();
    SparseMatrixXd SA2 = temp2.sparseView();
    VectorXd b2(6);
    b2 << 1, 1, 0, 0, 0, 0;

    res test2a = fr.reduce(SA2, b2, 6, true);
    REQUIRE((test2a.sparse_A.rows() == 4 && test2a.sparse_A.cols() == 5));
    res test2b = fr.reduce(SA2, b2, 6, false);
    REQUIRE((test2b.dense_A.rows() == 2 && test2b.dense_A.cols() == 1));

    MatrixXd A3(4,2);
    A3 << 1, 0, -1, 0, 0, 1, 0, -1;
    MatrixXd temp3 (4, 6);
    temp3 << A3, VectorXd::Ones(4).asDiagonal().toDenseMatrix();
    SparseMatrixXd SA3 = temp3.sparseView();
    VectorXd b3(4);
    b3 << 1, 0, 1, 0;

    res test3a = fr.reduce(SA3, b3, 4, true);
    REQUIRE((test3a.sparse_A.rows() == 4 && test3a.sparse_A.cols() == 6));
    res test3b = fr.reduce(SA3, b3, 4, false);
    REQUIRE((test3b.dense_A.rows() == 4 && test3b.dense_A.cols() == 2));

}

TEST_CASE( "Check Centering Algorithm", "[require]" ){
    SparseCenter sc;
    VectorXd simplex_x = sc.getInitialPoint(simplex.A, simplex.b, simplex.k);
    REQUIRE(simplex_x.rows() == 3);
    REQUIRE(simplex_x(0) == Catch::Approx(0.3333333).epsilon(0.01));
    REQUIRE(simplex_x(1) == Catch::Approx(0.3333333).epsilon(0.01));
    REQUIRE(simplex_x(2) == Catch::Approx(0.3333333).epsilon(0.01));

    VectorXd hc_x = sc.getInitialPoint(hc.A, hc.b, hc.k);
    REQUIRE(hc_x.rows() == 6);
    REQUIRE_THAT(hc_x(0), Catch::Matchers::WithinAbs(0, 0.0001));
    REQUIRE_THAT(hc_x(1), Catch::Matchers::WithinAbs(0, 0.0001));
    REQUIRE_THAT(hc_x(2), Catch::Matchers::WithinAbs(1, 0.0001));
    REQUIRE_THAT(hc_x(3), Catch::Matchers::WithinAbs(1, 0.0001));
    REQUIRE_THAT(hc_x(4), Catch::Matchers::WithinAbs(1, 0.0001));
    REQUIRE_THAT(hc_x(5), Catch::Matchers::WithinAbs(1, 0.0001));

    VectorXd birk_x = sc.getInitialPoint(birk.A, birk.b, birk.k);
    REQUIRE(birk_x.rows() == 4);
    REQUIRE_THAT(birk_x(0), Catch::Matchers::WithinAbs(0.5, 0.0001));
    REQUIRE_THAT(birk_x(1), Catch::Matchers::WithinAbs(0.5, 0.0001));
    REQUIRE_THAT(birk_x(2), Catch::Matchers::WithinAbs(0.5, 0.0001));
    REQUIRE_THAT(birk_x(3), Catch::Matchers::WithinAbs(0.5, 0.0001));

    DenseCenter dc; 

    MatrixXd A1 (4, 2);
    A1 << 1, 0, 0, 1, -1, 0, 0, -1;
    VectorXd b1 (4);
    b1 << 1, 1, 1, 1;

    VectorXd center1 = dc.getInitialPoint(A1, b1);
    REQUIRE_THAT(center1(0), Catch::Matchers::WithinAbs(0, 0.0001));
    REQUIRE_THAT(center1(1), Catch::Matchers::WithinAbs(0, 0.0001));

    MatrixXd A2 (3, 2);
    A2 << -1, 0, 0, -1, 1, 1;
    VectorXd b2 (3);
    b2 << 0, 0, 1;

    VectorXd center2 = dc.getInitialPoint(A2, b2);
    REQUIRE_THAT(center2(0), Catch::Matchers::WithinAbs(0.29, 0.01));
    REQUIRE_THAT(center2(1), Catch::Matchers::WithinAbs(0.29, 0.01));

}

TEST_CASE( "Check Weight Properties", "[require]" ){
    //Vaidya, John, DikinLS
    SparseVaidyaWalk vaidya_sparse(0.00001, 0.5);
    SparseDikinLSWalk dikinls_sparse(0.00001, 3, 0.001, 0.01, 10000);
    SparseJohnWalk john_sparse(0.00001, 0.5, 0.001, 0.01, 10000);

    VectorXd simplex_start (3);
    simplex_start << 0.33, 0.34, 0.33;
    SparseMatrixXd w = dikinls_sparse.generateWeight(simplex_start, simplex.A, simplex.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(2, 0.01));
    w = john_sparse.generateWeight(simplex_start, simplex.A, simplex.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(3, 0.01));
    w = vaidya_sparse.generateWeight(simplex_start, simplex.A, simplex.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(4, 0.01));

    VectorXd hc_start (6);
    hc_start << 0, 0, 1, 1, 1, 1;
    w = dikinls_sparse.generateWeight(hc_start, hc.A, hc.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(2, 0.01));
    w = john_sparse.generateWeight(hc_start, hc.A, hc.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(3, 0.01));
    w = vaidya_sparse.generateWeight(hc_start, hc.A, hc.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(4, 0.01));

    VectorXd birk_start (4);
    birk_start << 0.5, 0.5, 0.5, 0.5;
    w = dikinls_sparse.generateWeight(birk_start, birk.A, birk.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(1, 0.01));
    w = john_sparse.generateWeight(birk_start, birk.A, birk.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(1.5, 0.01));
    w = vaidya_sparse.generateWeight(birk_start, birk.A, birk.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(2, 0.01));

    FacialReduction fr;
    DenseCenter dc;
    res simplex_dense = fr.reduce(simplex.A, simplex.b, simplex.k, false);
    res hc_dense = fr.reduce(hc.A, hc.b, hc.k, false);
    res birk_dense = fr.reduce(birk.A, birk.b, birk.k, false);
    VectorXd sd_x = dc.getInitialPoint(simplex_dense.dense_A, simplex_dense.dense_b);
    VectorXd hc_x = dc.getInitialPoint(hc_dense.dense_A, hc_dense.dense_b);
    VectorXd birk_x = dc.getInitialPoint(birk_dense.dense_A, birk_dense.dense_b);

    JohnWalk john(0.5, 0.001, 0.01, 10000);
    DikinLSWalk dikinls(0.5, 0.001, 0.01, 10000);
    VaidyaWalk vaidya(0.5);

    john.generateWeight(sd_x, simplex_dense.dense_A, simplex_dense.dense_b);
    dikinls.generateWeight(sd_x, simplex_dense.dense_A, simplex_dense.dense_b);
    vaidya.generateWeight(sd_x, simplex_dense.dense_A, simplex_dense.dense_b);
    REQUIRE_THAT(dikinls.weights.diagonal().sum(), Catch::Matchers::WithinAbs(2, 0.01));
    REQUIRE_THAT(john.weights.diagonal().sum(), Catch::Matchers::WithinAbs(3, 0.01));
    REQUIRE_THAT(vaidya.weights.diagonal().sum(), Catch::Matchers::WithinAbs(4, 0.01));

    john.generateWeight(hc_x, hc_dense.dense_A, hc_dense.dense_b);
    dikinls.generateWeight(hc_x, hc_dense.dense_A, hc_dense.dense_b);
    vaidya.generateWeight(hc_x, hc_dense.dense_A, hc_dense.dense_b);
    REQUIRE_THAT(dikinls.weights.diagonal().sum(), Catch::Matchers::WithinAbs(2, 0.01));
    REQUIRE_THAT(john.weights.diagonal().sum(), Catch::Matchers::WithinAbs(3, 0.01));
    REQUIRE_THAT(vaidya.weights.diagonal().sum(), Catch::Matchers::WithinAbs(4, 0.01));

    john.generateWeight(birk_x, birk_dense.dense_A, birk_dense.dense_b);
    dikinls.generateWeight(birk_x, birk_dense.dense_A, birk_dense.dense_b);
    vaidya.generateWeight(birk_x, birk_dense.dense_A, birk_dense.dense_b);
    REQUIRE_THAT(dikinls.weights.diagonal().sum(), Catch::Matchers::WithinAbs(1, 0.01));
    REQUIRE_THAT(john.weights.diagonal().sum(), Catch::Matchers::WithinAbs(1.5, 0.01));
    REQUIRE_THAT(vaidya.weights.diagonal().sum(), Catch::Matchers::WithinAbs(2, 0.01));

}

TEST_CASE( "Test All Dense Combinations", "[require]" ){
    JohnWalk john(0.5, 0.001, 0.01, 100);
    DikinLSWalk dikinls(3.0, 0.001, 0.01, 100);
    VaidyaWalk vaidya(0.5);
    DikinWalk dikin(0.5);
    BallWalk ball(0.5);
    HitAndRun hitrun(0.001, 0.5);

    MatrixXd walk_res = denseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &john);
    walk_res = denseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &dikinls);
    walk_res = denseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &vaidya);
    walk_res = denseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &dikin);
    walk_res = denseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &ball);
    walk_res = denseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &hitrun);

    walk_res = denseFullWalkRun(hc.A, hc.b, hc.k, 100, &john);
    walk_res = denseFullWalkRun(hc.A, hc.b, hc.k, 100, &dikinls);
    walk_res = denseFullWalkRun(hc.A, hc.b, hc.k, 100, &vaidya);
    walk_res = denseFullWalkRun(hc.A, hc.b, hc.k, 100, &dikin);
    walk_res = denseFullWalkRun(hc.A, hc.b, hc.k, 100, &ball);
    walk_res = denseFullWalkRun(hc.A, hc.b, hc.k, 100, &hitrun);

    walk_res = denseFullWalkRun(birk.A, birk.b, birk.k, 100, &john);
    walk_res = denseFullWalkRun(birk.A, birk.b, birk.k, 100, &dikinls);
    walk_res = denseFullWalkRun(birk.A, birk.b, birk.k, 100, &vaidya);
    walk_res = denseFullWalkRun(birk.A, birk.b, birk.k, 100, &dikin);
    walk_res = denseFullWalkRun(birk.A, birk.b, birk.k, 100, &ball);
    walk_res = denseFullWalkRun(birk.A, birk.b, birk.k, 100, &hitrun);
}

TEST_CASE( "Test All Sparse Combinations", "[require]" ){
    SparseJohnWalk john(0.00001, 0.5, 0.001, 0.01, 100);
    SparseDikinLSWalk dikinls(0.00001, 3.0, 0.001, 0.01, 100);
    SparseVaidyaWalk vaidya(0.00001, 0.5);
    SparseDikinWalk dikin(0.00001, 0.5);
    SparseBallWalk ball(0.5);
    SparseHitAndRun hitrun(0.01, 0.5);

    MatrixXd walk_res = sparseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &john);
    walk_res = sparseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &dikinls);
    walk_res = sparseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &vaidya);
    walk_res = sparseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &dikin);
    walk_res = sparseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &ball);
    walk_res = sparseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &hitrun);

    walk_res = sparseFullWalkRun(hc.A, hc.b, hc.k, 100, &john);
    walk_res = sparseFullWalkRun(hc.A, hc.b, hc.k, 100, &dikinls);
    walk_res = sparseFullWalkRun(hc.A, hc.b, hc.k, 100, &vaidya);
    walk_res = sparseFullWalkRun(hc.A, hc.b, hc.k, 100, &dikin);
    walk_res = sparseFullWalkRun(hc.A, hc.b, hc.k, 100, &ball);
    walk_res = sparseFullWalkRun(hc.A, hc.b, hc.k, 100, &hitrun);

    walk_res = sparseFullWalkRun(birk.A, birk.b, birk.k, 100, &john);
    walk_res = sparseFullWalkRun(birk.A, birk.b, birk.k, 100, &dikinls);
    walk_res = sparseFullWalkRun(birk.A, birk.b, birk.k, 100, &vaidya);
    walk_res = sparseFullWalkRun(birk.A, birk.b, birk.k, 100, &dikin);
    walk_res = sparseFullWalkRun(birk.A, birk.b, birk.k, 100, &ball);
    walk_res = sparseFullWalkRun(birk.A, birk.b, birk.k, 100, &hitrun);
}