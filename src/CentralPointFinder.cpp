#include "CentralPointFinder.hpp"


VectorXd CentralPointFinder::getInitialPoint(MatrixXd A, VectorXd b){
    int n = A.rows();
    int d = A.cols();


    VectorXd y = VectorXd::Zero(d + 1);
    y(d) = ss_;

    VectorXd c = VectorXd::Zero(d + 1);
    c(d) = 1;

    VectorXd ones = VectorXd::Ones(n) * -1;

    for(int i = 0; i < n; i++){
        double norm = A.row(i).norm();
        A.row(i) = A.row(i)/norm;
        b(i) = b(i)/norm;
    }

    MatrixXd A_tilde(A.rows(), A.cols() + 1);
    A_tilde << A, ones;
    double t = ts_;

    VectorXd prev_grad = VectorXd::Zero(d + 1);

    while (t < te_){
        VectorXd slack = (b - (A_tilde * y)).cwiseInverse();

        VectorXd grad = c + A_tilde.transpose() * slack/t;
        MatrixXd slack_mat = slack.asDiagonal().toDenseMatrix();
        MatrixXd hess = A_tilde.transpose() * slack_mat * slack_mat * A_tilde/t;
        VectorXd sol = hess.colPivHouseholderQr().solve(grad);

        VectorXd new_y = y - 0.5 * (sol);
        if ((prev_grad - grad).norm()/(prev_grad.norm() + err_term_) < grad_lim_){
            t *= 2;
        }
        prev_grad = grad;
        y = new_y;
    }

    VectorXd ans = VectorXd(d);
    for(int i = 0; i < d; i++){
        ans(i) = y(i);
    }
    return ans;
}