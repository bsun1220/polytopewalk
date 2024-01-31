#include "FacialReduction.hpp"

MatrixXd FacialReduction::equalConversion(const MatrixXd& A){
    MatrixXd A_tilde (A.rows(), A.cols() + A.rows());
    A_tilde << A, MatrixXd::Identity(A.rows(), A.rows());
    return A_tilde;
}

z_result FacialReduction::findZ(const MatrixXd& newA, const VectorXd& b, int x_dim){
    int n = newA.rows();
    int d = newA.cols();
    z_result ans;
    ans.found_sol = false;
    if (x_dim + 2 > newA.rows()){
        return ans;
    }
    // setting up ineqA with extra column for delta
    MatrixXd ineqA (newA.cols() - x_dim, newA.rows() + 1);
    for(int i = x_dim; i < newA.cols(); i++){
        VectorXd temp_row (newA.rows() + 1);
        temp_row << newA.col(i), 1;
        ineqA.row(i - x_dim) = -1 * temp_row;
    }

    for(int ind = x_dim; ind < d; ind++){
         //setting up eqA
        MatrixXd eqA(x_dim + 2, newA.rows());

        for(int i = 0; i < x_dim; i++){
            eqA.row(i) = newA.col(i);
        }
        eqA.row(x_dim + 1) = b;
        eqA.row(x_dim) = newA.col(ind);

        VectorXd eqb = VectorXd::Zero(eqA.cols());
        eqb(x_dim) = 1;

        VectorXd init;
        FullPivLU <MatrixXd> lu(makeFullRank(eqA));
        //if eqA is not invertible skip the step
        if (!lu.isInvertible()) continue;
        //create init with last entry as delta
        init = lu.solve(eqb);

        double delta = ((-1 * newA).transpose() * init).maxCoeff();
        VectorXd temp (init.rows() + 1);
        temp << init, delta;
        init = temp;
        
        //add one column for delta
        MatrixXd temp_eqA (eqA.rows(), eqA.cols() + 1);
        temp_eqA << eqA, VectorXd::Zero(eqA.rows());
        eqA = temp_eqA;
        string name = "var_set1";
        if(init(init.rows() - 1) <= 0){
            init = init.head(init.rows() - 1);
            ans.found_sol = true; 
            ans.z = (newA.transpose() * init);
            return ans;
        }

        Problem lp;
        lp.AddVariableSet(make_shared<ExVariables>(n + 1, name, init));
        lp.AddConstraintSet(make_shared<ExConstraint1>(ineqA.rows(),name,ineqA));
        lp.AddConstraintSet(make_shared<ExConstraint2>(eqA.rows(),name,eqA, eqb));
        lp.AddCostSet(make_shared<ExCost>(name));
        
        IpoptSolver ipopt;
        ipopt.SetOption("print_level", 0);
        ipopt.SetOption("sb", "yes");
        ipopt.Solve(lp);

        VectorXd sol = lp.GetOptVariables()->GetValues();
        if (sol(sol.rows() - 1) <= 0){
            sol = sol.head(sol.rows() - 1);
            ans.found_sol = true; 
            ans.z = (newA.transpose() * sol);
            return ans;
        }
    }
    
    return ans;
}

MatrixXd FacialReduction::pickV(const VectorXd& z, int x_dim){
    int d = z.rows();
    vector<int> indices;
    for(int i = 0; i < x_dim; i++){
        indices.push_back(i); 
    }

    for(int i = x_dim; i < d; i++){
         if(z(i) <= 1.0E-8) indices.push_back(i);
    }

    MatrixXd matrix (indices.size(), d);
    for(int i = 0; i < indices.size(); i++){
        VectorXd row = VectorXd::Zero(d);
        row(indices[i]) = 1; 
        matrix.row(i) = row;
    }
    return matrix.transpose();
}

MatrixXd FacialReduction::pickP(const MatrixXd& AV){
    vector<int> lst;
    HouseholderQR <MatrixXd> qr(AV.cols(), AV.rows());
    qr.compute(AV.transpose());
    MatrixXd r = qr.matrixQR().triangularView<Eigen::Upper>();

    for(int i = 0; i < min(r.cols(), r.rows()); i++){
        if (abs(r.coeffRef(i, i)) > 1.0E-8){
            lst.push_back(i);
        }
    }

    MatrixXd proj (lst.size(), AV.rows());
    for(int i = 0; i < lst.size(); i++){
        VectorXd row = VectorXd::Zero(AV.rows());
        row(lst[i]) = 1; 
        proj.row(i) = row; 
    }
    return proj;
}

fr_result FacialReduction::entireFacialReductionStep(MatrixXd A, VectorXd b, int x_dim){
    z_result z_ans = findZ(A, b, x_dim);
    if(!z_ans.found_sol){
        fr_result f_ans;
        f_ans.A = A;
        f_ans.b = b;
        return f_ans;
    }
    MatrixXd V = pickV(z_ans.z, x_dim);
    MatrixXd AV = A * V;
    MatrixXd P = pickP(AV);

    A = P * AV;
    b = P * b;
    return entireFacialReductionStep(A, b, x_dim);
}

problem_result FacialReduction::reduce(MatrixXd A, VectorXd b){
    MatrixXd newA = equalConversion(A);
    int x_dim = A.cols();
    fr_result res = entireFacialReductionStep(newA, b, x_dim);
    if (res.A.cols() == newA.cols() && res.A.rows() == newA.rows()){
        problem_result ans;
        ans.reduced = false; 
        ans.reduced_A = A;
        ans.reduced_b = b;
        return ans;
    }

    HouseholderQR <MatrixXd> qr(res.A.cols(), res.A.rows());
    
    qr.compute(res.A.transpose());
    MatrixXd Q = qr.householderQ();
    MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();
    int d = R.rows();
    int n = R.cols();

    MatrixXd newR = R.block(0, 0, R.cols(), R.cols());
    VectorXd z1 = newR.transpose().inverse() * res.b;
    
    MatrixXd Q1 = Q.block(0, 0, Q.rows(), n);
    MatrixXd Q2 = Q.block(0, n, Q.rows(), d - n);

    MatrixXd reduced_A = -1 * Q2.block(x_dim, 0, d - x_dim, d - n);
    VectorXd reduced_b = (Q1 * z1).tail(d - x_dim);

    problem_result ans;
    ans.reduced = true; 
    ans.reduced_A = reduced_A;
    ans.reduced_b = reduced_b;
    
    ans.z1 = z1;
    ans.Q = Q;
    return ans;
    
}