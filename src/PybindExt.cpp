#include <pybind11/pybind11.h>
#include "PolytopeWalk.hpp"
#include <pybind11/eigen.h>
namespace py = pybind11;

MatrixXd generateWalk(MatrixXd& A, VectorXd& b, string walk_type, float r, int num_sim, float step_size, int max_iter, float grad_lim, double ss, double ts, double te, float err){
    
    CentralPointFinder cpf (ss, ts, te, err, grad_lim);
    FacialReduction fr;

    if (walk_type == "dikin"){
        DikinWalk d;
        return fullWalkRun(A, b, r, num_sim, d, fr, cpf);
    } else if (walk_type == "vaidya"){
        VaidyaWalk v;
        return fullWalkRun(A, b, r, num_sim, v, fr, cpf);
    } else if (walk_type == "dikinls"){
        DikinLSWalk dl(step_size, max_iter, grad_lim);
        return fullWalkRun(A, b, r, num_sim, dl, fr, cpf);
    } else if (walk_type == "john"){
        JohnWalk j(step_size, max_iter, grad_lim);
        return fullWalkRun(A, b, r, num_sim, j, fr, cpf);
    } else if (walk_type == "ball"){
        BallWalk ball;
        return fullWalkRun(A, b, r, num_sim, ball, fr, cpf);
    } else if (walk_type == "hitrun"){
        HitAndRunWalk hr (err);
        return fullWalkRun(A, b, r, num_sim, hr, fr, cpf);
    }
    return A;
}


PYBIND11_MODULE(polytopewalk, m) {
    m.doc() = "pybind11 poly walk example";
    m.def("generateWalk", &generateWalk, "central function", 
    py::arg("A"), py::arg("b"), py::arg("walk_type"), py::arg("r"), py::arg("num_sim"), 
    py::arg("step_size") = 0.1, py::arg("max_iter") = 100, py::arg("grad_lim") = 0.01, py::arg("ss") = 10000, 
    py::arg("ts") = 0.00001, py::arg("te") =10000, py::arg("err") = 0.01);

}
