#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <randomwalk/PolytopeWalk.hpp>
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

template <class RandomWalkBase = RandomWalk> class PyRandomWalk : public RandomWalkBase {
public:
    using RandomWalkBase::RandomWalkBase; // Inherit constructors
    MatrixXd generateCompleteWalk(const int num_steps, VectorXd& x, const MatrixXd& A, const VectorXd& b) override{
            PYBIND11_OVERRIDE_PURE(
                MatrixXd,
                RandomWalkBase,
                generateCompleteWalk,
                num_steps,
                x,
                A,
                b
            );
    }
};
template <class BarrierWalkBase = BarrierWalk> class PyBarrierWalk: public PyRandomWalk<BarrierWalkBase> {
public:
    using PyRandomWalk<BarrierWalkBase>::PyRandomWalk; // Inherit constructors
    // Override PyAnimal's pure virtual go() with a non-pure one:
    MatrixXd generateCompleteWalk(const int num_steps, VectorXd& x, const MatrixXd& A, const VectorXd& b) override{
            PYBIND11_OVERRIDE(
                MatrixXd,
                BarrierWalkBase,
                generateCompleteWalk,
                num_steps,
                x,
                A,
                b
            );
    }
    void generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b) override{
            PYBIND11_OVERRIDE(
                void,
                BarrierWalkBase,
                generateWeight,
                x,
                A,
                b
            );
        }
    void setDistTerm(int d, int n) override{
            PYBIND11_OVERRIDE(
                void,
                BarrierWalkBase,
                setDistTerm,
                d,
                n
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
    
    py::class_<RandomWalk, PyRandomWalk<>>(m, "RandomWalk")
        .def(py::init<>())
        .def("generateCompleteWalk", &RandomWalk::generateCompleteWalk);
    
    py::class_<BallWalk, RandomWalk>(m, "BallWalk")
        .def(py::init<const float>(), py::arg("r_p") = 0.3);
    
    py::class_<HitAndRunWalk, RandomWalk>(m, "HitAndRunWalk")
        .def(py::init<const float, const float>(), 
        py::arg("err_p") = 0.01, py::arg("r") = 0.1);

    py::class_<BarrierWalk, RandomWalk, PyBarrierWalk<>>(m, "BarrierWalk")
        .def(py::init<const float>(), py::arg("rp") = 1)
        .def("generateWeight", &BarrierWalk::generateWeight)
        .def("generateCompleteWalk", &RandomWalk::generateCompleteWalk);
    
    py::class_<DikinWalk, BarrierWalk, PyBarrierWalk<DikinWalk>>(m, "DikinWalk")
        .def(py::init<const float>(), py::arg("rp") = 1);
    
    py::class_<VaidyaWalk, BarrierWalk, PyBarrierWalk<VaidyaWalk>>(m, "VaidyaWalk")
        .def(py::init<const float>(), py::arg("rp") = 1);
    
    py::class_<DikinLSWalk, BarrierWalk, PyBarrierWalk<DikinLSWalk>>(m, "DikinLSWalk")
        .def(py::init<const float, const int, const float, const float>(), 
        py::arg("ss") = 0.1, py::arg("mi") = 10000, py::arg("gl") = 0.1, 
        py::arg("rp") = 1);
    
    py::class_<JohnWalk, BarrierWalk, PyBarrierWalk<JohnWalk>>(m, "JohnWalk")
        .def(py::init<const float, const int, const float, const float>(), 
        py::arg("ss") = 0.1, py::arg("mi") = 10000, py::arg("gl") = 0.1, 
        py::arg("rp") = 1);
    
    py::class_<problem_result>(m, "problem_result")
        .def_readwrite("reduced_A", &problem_result::reduced_A)
        .def_readwrite("reduced_b", &problem_result::reduced_b)
        .def_readwrite("reduced", &problem_result::reduced)
        .def_readwrite("z1", &problem_result::z1)
        .def_readwrite("Q", &problem_result::Q);

    py::class_<Reducer, PyReducer>(m, "Reducer")
        .def(py::init<>())
        .def("reduce", &Reducer::reduce);
    
    py::class_<FacialReduction, Reducer>(m, "FacialReduction")
        .def(py::init<>());

}