#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <constraintwalk/ConstraintWalk.hpp>
#include <constraintwalk/ConstraintBarrierWalk.hpp>
#include <constraintwalk/ConstraintDikinWalk.hpp>
#include <constraintwalk/ConstraintVaidyaWalk.hpp>
#include <constraintwalk/ConstraintHitRun.hpp>
#include <constraintwalk/ConstraintBallWalk.hpp>
#include <constraintwalk/SparseFacialReduction.hpp>
#include <constraintwalk/SparseCenter.hpp>
namespace py = pybind11;

template <class ConstraintWalkBase = ConstraintWalk> class PyConstraintWalk : public ConstraintWalkBase {
public:
    using ConstraintWalkBase::ConstraintWalkBase; // Inherit constructors
    MatrixXd generateCompleteWalk(
        const int num_steps, 
        const VectorXd& init, 
        const SparseMatrixXd& A, 
        const VectorXd& b, 
        int k) override
    {
            PYBIND11_OVERRIDE_PURE(
                MatrixXd,
                ConstraintWalkBase,
                generateCompleteWalk,
                num_steps,
                init,
                A,
                b,
                k
            );
    }
};
template <class ConsBarrierWalkBase = ConstraintBarrierWalk> class PyConsBarrierWalk: public PyConstraintWalk<ConsBarrierWalkBase> {
public:
    using PyConstraintWalk<ConsBarrierWalkBase>::PyConstraintWalk; 
    MatrixXd generateCompleteWalk(
        const int num_steps, 
        const VectorXd& init, 
        const SparseMatrixXd& A, 
        const VectorXd& b,
        int k) override
        {
            PYBIND11_OVERRIDE(
                MatrixXd,
                ConsBarrierWalkBase,
                generateCompleteWalk,
                num_steps,
                init,
                A,
                b,
                k
            );
    }

    SparseMatrixXd generateG(const VectorXd& x, const SparseMatrixXd& A, int k) override{
            PYBIND11_OVERRIDE(
                SparseMatrixXd,
                ConsBarrierWalkBase,
                generateG,
                x,
                A,
                k
            );
    }
    void setDistTerm(int d, int n) override{
            PYBIND11_OVERRIDE(
                void,
                ConsBarrierWalkBase,
                setDistTerm,
                d,
                n
            );
    }
};

PYBIND11_MODULE(conspolytopewalk, m) {
    m.doc() = "pybind11 sparse constraint walk library";
    
    py::class_<ConstraintWalk, PyConstraintWalk<>>(m, "ConstraintWalk")
        .def(py::init<const double>(), py::arg("err") = 0.000001)
        .def("generateCompleteWalk", &ConstraintWalk::generateCompleteWalk);
    
    py::class_<ConstraintBallWalk, ConstraintWalk>(m, "ConstraintBallWalk")
        .def(py::init<const double>(), py::arg("r") = 0.3);
    
    py::class_<ConstraintHitAndRun, ConstraintWalk>(m, "ConstraintHitAndRun")
        .def(py::init<const double, const double>(), 
        py::arg("err") = 0.01, py::arg("r") = 0.5);

    py::class_<ConstraintBarrierWalk, ConstraintWalk, PyConsBarrierWalk<>>(m, "ConstraintBarrierWalk")
        .def(py::init<const double, const double>(), py::arg("err") = 0.000001, py::arg("r") = 0.5)
        .def("generateG", &ConstraintBarrierWalk::generateG)
        .def("generateCompleteWalk", &ConstraintWalk::generateCompleteWalk);
    
    py::class_<ConstraintDikinWalk, ConstraintBarrierWalk, PyConsBarrierWalk<ConstraintDikinWalk>>(m, "ConstraintDikinWalk")
        .def(py::init<const double, const double>(), py::arg("err") = 0.000001, py::arg("r") = 0.5);
    
    py::class_<ConstraintVaidyaWalk, ConstraintBarrierWalk, PyConsBarrierWalk<ConstraintVaidyaWalk>>(m, "ConstraintVaidyaWalk")
        .def(py::init<const double, const double>(), py::arg("err") = 0.000001, py::arg("r") = 0.5);
    
    py::class_<SparseCenter>(m, "SparseCenter")
        .def(py::init<>())
        .def("getInitialPoint", &SparseCenter::getInitialPoint);
    
    py::class_<SparseFacialReduction>(m, "SparseFacialReduction")
        .def(py::init<>())
        .def("entireFacialReductionStep", &SparseFacialReduction::entireFacialReductionStep);
    
    py::class_<fr_res>(m, "fr_res")
        .def_readwrite("A", &fr_res::A)
        .def_readwrite("b", &fr_res::b);
}