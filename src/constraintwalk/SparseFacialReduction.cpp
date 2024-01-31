#include "SparseFacialReduction.hpp"

z_res SparseFacialReduction::findZ(const SparseMatrixXd& A, const VectorXd& b, int x_dim){
    int n = A.rows();
    int d = A.cols();

    z_res ans;
    ans.found_sol = false;

    if (x_dim + 2 > n){
        return ans;
    }

    SparseMatrixXd ineqA(n + 1, d-x_dim);
    for(int i = x_dim; i < d; i++){
        ineqA.col(i - x_dim) = -1 * A.col(i);
    }
    ineqA = ineqA.transpose();
    ineqA.col(n) = (-1 * VectorXd::Ones(d-x_dim)).sparseView();

    SparseMatrixXd eqA(n + 1, x_dim + 2);
    for(int j = 0; j < x_dim; j++){
        eqA.col(j) = A.col(j);
    }
    eqA.col(x_dim + 1) = b.sparseView();
    VectorXd eqb = VectorXd::Zero(eqA.cols());
    eqb(x_dim) = 1;
    
    for(int i = x_dim; i < d; i++){
        eqA.col(x_dim) = A.col(i);
        eqA = eqA.transpose();

        SparseQR<SparseMatrixXd, COLAMDOrdering<SparseMatrix<double>::StorageIndex>> solver (eqA.block(0, 0, x_dim + 2, n));
        VectorXd init = solver.solve(eqb);
        double delta = ((-1 * A).transpose() * init).maxCoeff();
        VectorXd temp (init.rows() + 1);
        temp << init, delta;
        init = temp; 

        if(!eqb.isApprox(eqA * init)){
            eqA = eqA.transpose();
            continue;
        }

        if(init(init.rows() - 1) <= 0){
            init = init.head(init.rows() - 1);
            ans.found_sol = true; 
            ans.z = (A.transpose() * init);
            return ans;
        }

        string name = "var_set1";
        Problem lp;
        lp.AddVariableSet(make_shared<SparseExVariables>(n + 1, name, init));
        lp.AddConstraintSet(make_shared<SparseExConstraint1>(ineqA.rows(),name,ineqA));
        lp.AddConstraintSet(make_shared<SparseExConstraint2>(eqA.rows(),name,eqA, eqb));
        lp.AddCostSet(make_shared<SparseExCost>(name));
        IpoptSolver ipopt;
        ipopt.SetOption("print_level", 0);
        ipopt.SetOption("sb", "yes");
        ipopt.Solve(lp);

        VectorXd sol = lp.GetOptVariables()->GetValues();
        if (sol(sol.rows() - 1) <= 0){
            sol = sol.head(sol.rows() - 1);
            ans.found_sol = true; 
            ans.z = (A.transpose() * sol);
            return ans;
        }
        eqA = eqA.transpose();

    }
    return ans;
}

SparseMatrixXd SparseFacialReduction::pickV(const VectorXd& z, int x_dim){
    int d = z.rows();
    
    vector<T> indices;
    for(int i = 0; i < x_dim; i++){
        indices.push_back(T(indices.size(), i, 1)); 
    }
    for(int i = x_dim; i < d; i++){
         if(z(i) <= 1.0E-8) indices.push_back(T(indices.size(), i, 1)); 
    }
    SparseMatrixXd mat(indices.size(), d);
    mat.setFromTriplets(indices.begin(), indices.end());
    return mat.transpose();
}

SparseMatrixXd SparseFacialReduction::pickP(const SparseMatrixXd& AV){

    SparseQR<SparseMatrixXd, NaturalOrdering<SparseMatrix<double>::StorageIndex>> solver;
    solver.compute(AV.transpose());
    SparseMatrixXd R = solver.matrixR();

    vector<T> indices;
    for (int i = 0; i < min(R.cols(), R.rows()); i++){
        if (abs(R.coeffRef(i, i)) > 1.0E-8){
            indices.push_back(T(indices.size(), solver.colsPermutation().indices()(i), 1));
        }
    }
    SparseMatrixXd proj (indices.size(), AV.rows());
    proj.setFromTriplets(indices.begin(), indices.end());
    return proj; 
}

fr_res SparseFacialReduction::entireFacialReductionStep(SparseMatrixXd A, VectorXd b, int x_dim){
    z_res z_ans = findZ(A, b, x_dim);
    if(!z_ans.found_sol){
        fr_res ans;
        ans.A = A;
        ans.b = b; 
        return ans; 
    }
    SparseMatrixXd V = pickV(z_ans.z, x_dim);
    SparseMatrixXd AV = A * V;
    SparseMatrixXd P = pickP(AV);
    A = P * AV;
    b = P * b; 
    return entireFacialReductionStep(A, b, x_dim);
}