#include <pybind11/pybind11.h>
#include "PolytopeWalk.hpp"
#include <pybind11/eigen.h>
namespace py = pybind11;

MatrixXd generateWalk(MatrixXd& A, VectorXd& b, string walk_type, float r, int num_sim, float step_size, int max_iter, float grad_lim, double ss, double ts, double te, float err){
    
    CentralPointFinder cpf (ss, ts, te, err, grad_lim);
    FacialReduction fr;

    if (walk_type == "dikin"){
        DikinWalk d(r);
        return fullWalkRun(A, b, num_sim, d, fr, cpf);
    } else if (walk_type == "vaidya"){
        VaidyaWalk v(r);
        return fullWalkRun(A, b, num_sim, v, fr, cpf);
    } else if (walk_type == "dikinls"){
        DikinLSWalk dl(step_size, max_iter, grad_lim, r);
        return fullWalkRun(A, b, num_sim, dl, fr, cpf);
    } else if (walk_type == "john"){
        JohnWalk j(step_size, max_iter, grad_lim, r);
        return fullWalkRun(A, b, num_sim, j, fr, cpf);
    } else if (walk_type == "ball"){
        BallWalk ball (r);
        return fullWalkRun(A, b, num_sim, ball, fr, cpf);
    } else if (walk_type == "hitrun"){
        HitAndRunWalk hr (err, r);
        return fullWalkRun(A, b, num_sim, hr, fr, cpf);
    }
    return A;
};


class PyInitializer : public Initializer{
    public:
        using Initializer::Initializer;

        VectorXd getInitialPoint(MatrixXd A, VectorXd b) override{
            PYBIND11_OVERRIDE_PURE(
                Eigen::VectorXd,
                Initializer,
                getInitialPoint,
                A, 
                b
            );
        }
};

class PyRandomWalk : public RandomWalk{
    public:
        using RandomWalk::RandomWalk;

        MatrixXd generateCompleteWalk(const int num_steps, VectorXd& x, const MatrixXd& A, const VectorXd& b) override{
            PYBIND11_OVERRIDE_PURE(
                MatrixXd,
                RandomWalk,
                generateCompleteWalk,
                num_steps,
                x,
                A,
                b
            );
        }
};

class PyBarrierWalk : public BarrierWalk{
    public:
        using BarrierWalk::BarrierWalk;
        void generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b) override{
            PYBIND11_OVERRIDE_PURE(
                void,
                BarrierWalk,
                generateWeight,
                x,
                A,
                b
            );
        }

};


PYBIND11_MODULE(polytopewalk, m) {
    m.doc() = "pybind11 polytopwalk library";
    m.def("generateWalk", &generateWalk, "central function", 
    py::arg("A"), py::arg("b"), py::arg("walk_type"), py::arg("r"), py::arg("num_sim"), 
    py::arg("step_size") = 0.1, py::arg("max_iter") = 100, py::arg("grad_lim") = 0.01, py::arg("ss") = 10000, 
    py::arg("ts") = 0.00001, py::arg("te") =10000, py::arg("err") = 0.01);
    
    py::class_<Initializer, PyInitializer>(m, "Initializer")
        .def(py::init<>())
        .def("getInitialPoint", &Initializer::getInitialPoint);
    py::class_<CentralPointFinder, Initializer>(m, "CentralPointFinder")
        .def(py::init<const double, const double, 
            const double, const double, const double >(),
            py::arg("ss") = 10000, py::arg("s") = 0.00001, 
            py::arg("te") = 10000, py::arg("err_term") = 0.000001,
            py::arg("grad_lim") = 0.01);
    
    py::class_<RandomWalk, PyRandomWalk>(m, "RandomWalk")
        .def(py::init<>())
        .def("generateCompleteWalk", &RandomWalk::generateCompleteWalk);
    
    py::class_<BallWalk, RandomWalk>(m, "BallWalk")
        .def(py::init<const float>(), py::arg("r_p") = 0.3);
    
    py::class_<HitAndRunWalk, RandomWalk>(m, "HitAndRunWalk")
        .def(py::init<const float, const float>(), 
        py::arg("err_p") = 0.01, py::arg("r") = 0.1);

    py::class_<BarrierWalk, PyBarrierWalk>(m, "BarrierWalk")
        .def(py::init<const float>(), py::arg("rp") = 1)
        .def("generateWeight", &BarrierWalk::generateWeight)
        .def("generateCompleteWalk", &RandomWalk::generateCompleteWalk);
    
    py::class_<DikinWalk, BarrierWalk>(m, "DikinWalk")
        .def(py::init<const float>(), py::arg("rp") = 1);
    
    py::class_<VaidyaWalk, BarrierWalk>(m, "VaidyaWalk")
        .def(py::init<const float>(), py::arg("rp") = 1);
    
    py::class_<DikinLSWalk, BarrierWalk>(m, "DikinLSWalk")
        .def(py::init<const float, const int, const float, const float>(), 
        py::arg("ss") = 0.1, py::arg("mi") = 100, py::arg("mi") = 0.1, 
        py::arg("rp") = 1);
    
    py::class_<JohnWalk, BarrierWalk>(m, "JohnWalk")
        .def(py::init<const float, const int, const float, const float>(), 
        py::arg("ss") = 0.1, py::arg("mi") = 100, py::arg("mi") = 0.1, 
        py::arg("rp") = 1);

}