#include "FacialReductionF2.hpp"

MatrixXd FacialReductionF2::equalConversion(const MatrixXd& A){
    MatrixXd A_tilde (A.rows(), A.cols() + A.rows());
    A_tilde << A, MatrixXd::Identity(A.rows(), A.rows());
    return A_tilde;
}

z_result FacialReductionF2::findZ(const MatrixXd& newA, const VectorXd& b, int x_dim){
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
        eqA = makeFullRank(eqA);
        FullPivLU <MatrixXd> lu(eqA);
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
            VectorXd ans_v (sol.rows() - 1);
            for(int i = 0; i < ans_v.rows(); i++){
                ans_v(i) = sol(i);
            }
            ans.found_sol = true; 
            ans.z = (newA.transpose() * ans_v);
            return ans;
        }
    }
    
    return ans;
}

MatrixXd FacialReductionF2::pickV(const VectorXd& z, int x_dim){
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

MatrixXd FacialReductionF2::pickP(const MatrixXd& AV){
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

fr_result FacialReductionF2::entireFacialReductionStep(MatrixXd A, VectorXd b, int x_dim){
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

fr_result FacialReductionF2::reduceSampling(const MatrixXd& M, const VectorXd& b, int delta_dim){
    MatrixXd M_inv = M.inverse();
    int row = M_inv.rows() - delta_dim;
    int gamma = M_inv.rows() - b.rows();

    MatrixXd temp1(M_inv.rows() - row, M_inv.cols() - gamma);
    for(int i = 0; i < temp1.rows(); i++){
        for(int k = 0; k < temp1.cols(); k++){
            temp1(i, k) = M_inv(i + row, k);
        }
    }
    MatrixXd temp2(M_inv.rows() - row, gamma);
    for(int i = 0; i < temp2.rows(); i++){
        for(int k = 0; k < temp2.cols(); k++){
            temp2(i, k) = M_inv(i + row, k + M_inv.cols() - gamma);
        }
    }
    fr_result ans;
    ans.A = -1 * temp2;
    ans.b = temp1 * b;
    return ans;
}

problem_result FacialReductionF2::reduce(MatrixXd A, VectorXd b){
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

    MatrixXd M = makeFullRank(res.A);
    
    fr_result reduced = reduceSampling(M, res.b, res.A.cols() - x_dim);

    cout << "reduced_A" << endl;
    cout << reduced.A << endl;
    cout << "reduced_b" << endl;
    cout << reduced.b << endl;

    HouseholderQR <MatrixXd> qr(res.A.cols(), res.A.rows());
    qr.compute(res.A.transpose());
    MatrixXd q = qr.householderQ();
    MatrixXd r = qr.matrixQR().triangularView<Eigen::Upper>();
    
    int d = r.rows();
    int n = r.cols();

    MatrixXd new_r (r.cols(), r.cols());
    for(int i = 0; i < r.cols(); i++){
        new_r.row(i) = r.row(i);
    }
    VectorXd z1 = new_r.transpose().inverse() * res.b;
    
    MatrixXd Q1 (d, n);
    MatrixXd Q2 (d, d-n);
    for(int i = 0; i < n; i++){
        Q1.col(i) = q.col(i);
    }
    for(int i = d; i < d - n; i++){
        Q2.col(i - d) = q.col(i);
    }

    problem_result ans;
    ans.reduced = true; 
    ans.reduced_A = reduced.A;
    ans.reduced_b = reduced.b;
    ans.b_tilde = res.b;
    ans.M = M;
    return ans;
    
}