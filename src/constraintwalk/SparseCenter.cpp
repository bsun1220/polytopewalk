#include "SparseCenter.hpp"

VectorXd SparseCenter::getInitialPoint(SparseMatrixXd A, VectorXd b, int k){
    SparseQR<SparseMatrixXd, COLAMDOrdering<SparseMatrix<double>::StorageIndex>> solver (A);
    VectorXd solved = solver.solve(b);
    double delta = -1 * solved.minCoeff();
    VectorXd init (A.cols() + 1);
    init << solved, delta; 
    SparseMatrixXd eqA(A.rows(), A.cols() + 1);
    for(int i = 0; i < A.cols(); i++){
        eqA.col(i) = A.col(i);
    }
    SparseMatrixXd ineqA(k, A.cols() + 1);
    vector<T> indices;
    for(int i = A.cols() - k; i < A.cols(); i++){
        indices.push_back(T(i, i, -1));
        indices.push_back(T(i, A.cols(), -1));
    }
    ineqA.setFromTriplets(indices.begin(), indices.end());

    int n = A.cols();
    string name = "var_set1";
    Problem lp;
    lp.AddVariableSet(make_shared<SparseExVariables>(n + 1, name, init));
    lp.AddConstraintSet(make_shared<SparseExConstraint1>(ineqA.rows(),name,ineqA));
    lp.AddConstraintSet(make_shared<SparseExConstraint2>(eqA.rows(),name,eqA, b));
    lp.AddCostSet(make_shared<SparseExCost>(name));
    IpoptSolver ipopt;
    ipopt.SetOption("print_level", 0);
    ipopt.SetOption("sb", "yes");
    ipopt.Solve(lp);
    VectorXd sol = lp.GetOptVariables()->GetValues();
    if (sol(sol.rows() - 1) <= 0){
        sol = sol.head(sol.rows() - 1);
        return sol; 
    }
    cout << "NO CENTER FOUND" << endl;
    return sol;
    
}