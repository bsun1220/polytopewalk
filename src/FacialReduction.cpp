#include "FacialReduction.hpp"

class ExVariables : public VariableSet{
    public:
    VectorXd x;

    ExVariables(int num_dim, string name, VectorXd init) : VariableSet(num_dim, name){
        x = init;
    } 
    void SetVariables(const VectorXd& val) override
    {
        x = val;
    }

    VectorXd GetValues() const override{
        return x;
    }

    VecBound GetBounds() const override{
        VecBound bounds(GetRows());
        for(int i = 0; i < x.rows(); i++){
            bounds.at(i) = NoBound;
        }
        return bounds;
    }

};

class ExConstraint1 : public ConstraintSet{
    public:
    MatrixXd A;
    ExConstraint1(int num_dim, string name, MatrixXd A_param) : ConstraintSet(num_dim, name){
        // A is d by n matrix
        A = A_param;
        
    }

    VectorXd GetValues() const override{
        VectorXd g(GetRows());
        VectorXd input = GetVariables()->GetComponent("var_set1")->GetValues();
        for(int i = 0; i < A.rows(); i++){
            VectorXd A_row = A.row(i);
            //A row is a n dimensional vector
            g(i) = A_row.transpose() * input; 
        }
        return g;
    }

    VecBound GetBounds() const override{
        VecBound b(GetRows());
        for(int i = 0; i < A.rows(); i++){
            b.at(i) = Bounds(-inf, 0);
        }
        return b;
    }
    void FillJacobianBlock (string var_set, Jacobian& jac_block) const override{}

};

class ExConstraint2 : public ConstraintSet{
    public:
    MatrixXd A;
    VectorXd b;
    ExConstraint2(int num_dim, string name, MatrixXd A_param, VectorXd b_param) : ConstraintSet(num_dim, name){
        A = A_param;
        b = b_param;

    }

    VectorXd GetValues() const override{
        VectorXd g(GetRows());
        VectorXd input = GetVariables()->GetComponent("var_set1")->GetValues();
        for(int i = 0; i < A.rows(); i++){
            VectorXd A_row = A.row(i);
            //A row is a n dimensional vector
            g(i) = A_row.transpose() * input;

            const double multiplier = std::pow(10.0, 10);
             g(i) = round(g(i) * multiplier) / multiplier;
        }
        return g;
    }

    VecBound GetBounds() const override{
        VecBound bound(GetRows());
        for(int i = 0; i < b.rows(); i++){
            bound.at(i) = Bounds(b(i),b(i));
        }
        return bound;
    }

    void FillJacobianBlock (string var_set, Jacobian& jac_block) const override{}

};

class ExCost : public CostTerm{
    public:
    ExCost(string name) : CostTerm(name) {
    }
    
    double GetCost() const override
    {
        VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();
        return x(x.rows() - 1);
    };
    void FillJacobianBlock (string var_set, Jacobian& jac_block) const override{}
};


MatrixXd FacialReduction::equalConversion(MatrixXd A){
    MatrixXd A_tilde (A.rows(), A.cols() + A.rows());
    A_tilde << A, MatrixXd::Identity(A.rows(), A.rows());
    return A_tilde;
}

z_result FacialReduction::findZ(MatrixXd newA, VectorXd b, int x_dim){
    int n = newA.rows();
    int d = newA.cols();

    // setting up ineqA with extra column for delta
    MatrixXd ineqA (newA.rows() + 1, newA.cols());
    for(int i = 0; i < newA.rows(); i++){
        ineqA.row(i) = -1 * newA.row(i);
    }
    VectorXd temp_row(newA.cols());
    for(int i = 0; i < newA.cols(); i++){
        temp_row(i) = -1;
    }

    ineqA.row(newA.rows()) = temp_row;
    MatrixXd temp_m = ineqA.transpose();
    ineqA = temp_m;

    for(int ind = x_dim; ind < d; ind++){
        //setting up eqA
        MatrixXd eqA(x_dim + 2, newA.rows());
        for(int i = 0; i < x_dim; i++){
                eqA.row(i) = newA.col(i);
        }
        eqA.row(x_dim) = newA.col(ind);
        eqA.row(x_dim + 1) = b;
        int b_length = eqA.cols() > x_dim + 2 ? eqA.cols() : x_dim + 2;
        VectorXd eqb = VectorXd::Zero(b_length);
        eqb(x_dim) = 1;

        VectorXd init;
        if (eqA.rows() <= eqA.cols()){
            eqA = makeFullRank(eqA);
            FullPivLU <MatrixXd> lu(eqA);
            //if eqA is not invertible skip the step
            if (!lu.isInvertible()) continue;
            //create init with last entry as delta
            init = eqA.inverse() * eqb;
        } else {
            init = eqA.colPivHouseholderQr().solve(eqb);
        }
        
        double delta = ((-1 * newA).transpose() * init).maxCoeff();
        VectorXd temp (init.rows() + 1);
        for(int i =0 ; i < init.rows(); i++){
            temp(i) = init(i);
        }
        temp(init.rows()) = delta;
        init = temp;

        //add one column for delta
        MatrixXd temp_eqA (eqA.rows(), eqA.cols() + 1);
        for(int i = 0; i < eqA.cols(); i++){
            temp_eqA.col(i) = eqA.col(i);
        }
       
        temp_eqA.col(eqA.cols()) = VectorXd::Zero(eqA.rows());
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
            z_result ans;
            ans.found_sol = true; 
            ans.z = (newA.transpose() * ans_v);
            return ans;
        }
    }
    z_result ans;
    ans.found_sol = false;
    return ans;
}

MatrixXd FacialReduction::facialReduction(VectorXd z){
    int d = z.rows();
    vector<int> indices;
    const double multiplier = std::pow(10.0, 5);
    for(int i = 0; i < d; i++){
        z(i) = round(z(i) * multiplier) / multiplier;
        if(z(i) <= 0) indices.push_back(i);
    }
    MatrixXd matrix (indices.size(), d);
    for(int i = 0; i < indices.size(); i++){
        VectorXd row = VectorXd::Zero(d);
        row(indices[i]) = 1; 
        matrix.row(i) = row;
    }
    return matrix.transpose();
}

fr_result FacialReduction::entireFacialReductionStep(MatrixXd A, VectorXd b, int x_dim){
    z_result z_ans = findZ(A, b, x_dim);
    if(!z_ans.found_sol){
        fr_result f_ans;
        f_ans.A = A;
        f_ans.b = b;
        return f_ans;
    }

    MatrixXd V = facialReduction(z_ans.z);
    MatrixXd AV = A * V;

    FullPivLU<MatrixXd> lu_decomp(AV.transpose());
    MatrixXd decomp = lu_decomp.image(AV.transpose());

    vector<int> lst; 
    for(int i = 0; i < AV.rows(); i++){
        for(int j = 0; j< decomp.cols(); j++){
            if(AV.row(i).transpose().isApprox(decomp.col(j))){
                lst.push_back(i);
            }
        }
    }

    MatrixXd proj (lst.size(), A.rows());
    for(int i = 0; i < lst.size(); i++){
        VectorXd row = VectorXd::Zero(A.rows());
        row(lst[i]) = 1; 
        proj.row(i) = row; 
    }
    A = proj * AV;
    b = proj * b;

    return entireFacialReductionStep(A, b, x_dim);
}

fr_result FacialReduction::reduceSampling(MatrixXd M, VectorXd b, int delta_dim){
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
    MatrixXd M = makeFullRank(res.A);
    fr_result reduced = reduceSampling(M, res.b, res.A.cols() - x_dim);

    problem_result ans;
    ans.reduced = true; 
    ans.reduced_A = reduced.A;
    ans.reduced_b = reduced.b;
    ans.b_tilde = res.b;
    ans.M = M;
    return ans;
    
}