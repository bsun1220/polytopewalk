#include "constraintwalk/SparseFacialReduction.hpp"
#include "utils/FacialReduction.hpp"
#include "constraintwalk/SparseCenter.hpp"

int main(){
    FacialReduction fr; 
    SparseFacialReduction sfr; 

    MatrixXd A1 (6, 3);
    A1 << 1, 1, 0, -1, -1, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;

    VectorXd b1(6);
    b1 << 1, -1, 1, 1, 1, 1;

    A1 = fr.equalConversion(A1);
    fr_result res = fr.entireFacialReductionStep(A1, b1, 3);

    SparseMatrixXd SA1 = A1.sparseView();
    fr_res sfr_res = sfr.entireFacialReductionStep(SA1, b1, 3);
    
    MatrixXd A2(6,3);
    A2 << 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;
    cout << "START" << endl;
    VectorXd b2(6);
    b2 << 1, 1, 0, 0, 0, 0;
    A2 = fr.equalConversion(A2); 
    res = fr.entireFacialReductionStep(A2, b2, 3);
    cout << "A" << endl;
    cout << res.A << endl;
    cout << "b" << endl;
    cout << res.b << endl;
    cout << "---" << endl;

    SparseMatrixXd SA2 = A2.sparseView();
    sfr_res = sfr.entireFacialReductionStep(SA2, b2, 3);
    cout << "A" << endl;
    cout << MatrixXd(sfr_res.A) << endl;
    cout << "b" << endl;
    cout << sfr_res.b << endl;


    MatrixXd prop_A (1,4);
    prop_A << 1,1,1,1;
    VectorXd b(1);
    b << 1; 
    SparseMatrixXd A = prop_A.sparseView();
    SparseCenter sp;
    VectorXd sol = sp.getInitialPoint(A, b, 4);
    cout << sol << endl;

}
