#include <pybind11/pybind11.h>
#include "PolytopeWalk.hpp"
#include <pybind11/eigen.h>
namespace py = pybind11;

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

class PyReducer : public Reducer{
    public:
        using Reducer::Reducer;
        problem_result reduce(MatrixXd A, VectorXd b) override{
            PYBIND11_OVERRIDE_PURE(
                problem_result,
                Reducer,
                reduce,
                A,
                b
            );
        }
};


PYBIND11_MODULE(polytopewalk, m) {
    m.doc() = "pybind11 polytopewalk library";

    
    m.def("fullWalkRun", &fullWalkRun, "Central Function", py::arg("A"), 
    py::arg("b"), py::arg("num_sim"), py::arg("walk"), py::arg("reducer"), 
    py::arg("initializer"));
    
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
    
    py::class_<problem_result>(m, "problem_result")
        .def_readwrite("reduced_A", &problem_result::reduced_A)
        .def_readwrite("reduced_b", &problem_result::reduced_b)
        .def_readwrite("reduced", &problem_result::reduced)
        .def_readwrite("b_tilde", &problem_result::b_tilde)
        .def_readwrite("M", &problem_result::M);

    py::class_<Reducer, PyReducer>(m, "Reducer")
        .def(py::init<>())
        .def("reduce", &Reducer::reduce);
    
    py::class_<FacialReduction, Reducer>(m, "FacialReduction")
        .def(py::init<>());

}