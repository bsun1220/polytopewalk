#include "HitRunWalk.hpp"

double HitAndRunWalk::distance(VectorXd& x, VectorXd&y){
    return (x - y).norm();
}

double HitAndRunWalk::binarySearch(VectorXd direction, VectorXd x){

    VectorXd farth = x + r * direction;
    double dist = 0; 

    while(true){
        dist = distance(x, farth);
        farth = x + 2 * dist * direction; 
        if (!acceptReject(farth, A, b)){
            break; 
        }
    }
    VectorXd left = x;
    VectorXd right = farth;
    VectorXd mid = (x + farth)/2;

    while (distance(left, right) > err || ! acceptReject(mid, A, b)){
        mid = (left + right)/2; 
        if (acceptReject(mid, A, b)){
            left = mid; 
        } else {
            right = mid; 
        }
    }
    return distance(mid, x);
}

MatrixXd HitAndRunWalk::generateCompleteWalk(int num_steps, VectorXd x){
    int n = x.rows(); 
    MatrixXd results = MatrixXd::Zero(num_steps, n);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);
    for (int i = 0; i < num_steps; i++){
        VectorXd new_direct = generateGaussianRVNorm(n);
        double pos_side = binarySearch(new_direct, x);
        double neg_side = binarySearch(new_direct * -1, x) * -1;
        float val = dis(gen);
        double random_point = val * (pos_side - neg_side) + neg_side; 
        x = random_point * new_direct + x; 
        results.row(i) = x; 
    }
    return results;
}