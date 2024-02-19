#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "utils/FullWalkRun.hpp"
namespace py = pybind11;

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


template <class SparseRandomWalkBase = SparseRandomWalk> class PySparseRandomWalk : public SparseRandomWalkBase {
public:
    using SparseRandomWalkBase::SparseRandomWalkBase; // Inherit constructors
    MatrixXd generateCompleteWalk(
        const int num_steps, 
        const VectorXd& init, 
        const SparseMatrixXd& A, 
        const VectorXd& b, 
        int k) override
    {
            PYBIND11_OVERRIDE_PURE(
                MatrixXd,
                SparseRandomWalkBase,
                generateCompleteWalk,
                num_steps,
                init,
                A,
                b,
                k
            );
    }
};
template <class SparseBarrierWalkBase = SparseBarrierWalk> class PySparseBarrierWalk: public PySparseRandomWalk<SparseBarrierWalkBase> {
public:
    using PySparseRandomWalk<SparseBarrierWalkBase>::PySparseRandomWalk; 
    MatrixXd generateCompleteWalk(
        const int num_steps, 
        const VectorXd& init, 
        const SparseMatrixXd& A, 
        const VectorXd& b,
        int k) override
        {
            PYBIND11_OVERRIDE(
                MatrixXd,
                SparseBarrierWalkBase,
                generateCompleteWalk,
                num_steps,
                init,
                A,
                b,
                k
            );
    }

    SparseMatrixXd generateWeight(const VectorXd& x, const SparseMatrixXd& A, int k) override{
            PYBIND11_OVERRIDE(
                SparseMatrixXd,
                SparseBarrierWalkBase,
                generateWeight,
                x,
                A,
                k
            );
    }
    void setDistTerm(int d, int n) override{
            PYBIND11_OVERRIDE(
                void,
                SparseBarrierWalkBase,
                setDistTerm,
                d,
                n
            );
    }
};


PYBIND11_MODULE(polytopewalk, m) {
    m.doc() = "pybind11 polytopewalk library";

    m.def("denseFullWalkRun", &denseFullWalkRun, "Dense Central Function", py::arg("A"), 
    py::arg("b"), py::arg("k"), py::arg("num_sim"), py::arg("walk"));

    m.def("sparseFullWalkRun", &sparseFullWalkRun, "Sparse Central Function", py::arg("A"), 
    py::arg("b"), py::arg("k"), py::arg("num_sim"), py::arg("walk"));

    auto m_dense = m.def_submodule("dense", "Dense Module");
    auto m_sparse = m.def_submodule("sparse", "Sparse Module");

    py::class_<DenseCenter>(m_dense, "DenseCenter")
        .def(py::init<>())
        .def("getInitialPoint", &DenseCenter::getInitialPoint);
    
    py::class_<SparseCenter>(m_sparse, "SparseCenter")
        .def(py::init<>())
        .def("getInitialPoint", &SparseCenter::getInitialPoint);
    
    py::class_<RandomWalk, PyRandomWalk<>>(m_dense, "RandomWalk")
        .def(py::init<>())
        .def("generateCompleteWalk", &RandomWalk::generateCompleteWalk);
    
    py::class_<BallWalk, RandomWalk>(m_dense, "BallWalk")
        .def(py::init<const double>(), py::arg("r") = 0.3);
    
    py::class_<HitAndRun, RandomWalk>(m_dense, "HitAndRun")
        .def(py::init<const double, const double>(), 
        py::arg("err") = 0.01, py::arg("r") = 0.1);

    py::class_<BarrierWalk, RandomWalk, PyBarrierWalk<>>(m_dense, "BarrierWalk")
        .def(py::init<const double>(), py::arg("r") = 1)
        .def("generateWeight", &BarrierWalk::generateWeight)
        .def("generateCompleteWalk", &RandomWalk::generateCompleteWalk);
    
    py::class_<DikinWalk, BarrierWalk, PyBarrierWalk<DikinWalk>>(m_dense, "DikinWalk")
        .def(py::init<const double>(), py::arg("r") = 1);
    
    py::class_<VaidyaWalk, BarrierWalk, PyBarrierWalk<VaidyaWalk>>(m_dense, "VaidyaWalk")
        .def(py::init<const double>(), py::arg("r") = 1);
    
    py::class_<DikinLSWalk, BarrierWalk, PyBarrierWalk<DikinLSWalk>>(m_dense, "DikinLSWalk")
        .def(py::init<const double, const double, const double, const int>(), 
        py::arg("r") = 1.0, py::arg("g_lim") = 0.05, py::arg("step_size") = 1.0, 
        py::arg("max_iter") = 1000);
    
    py::class_<JohnWalk, BarrierWalk, PyBarrierWalk<JohnWalk>>(m_dense, "JohnWalk")
        .def(py::init<const double, const double, const double, const int>(), 
        py::arg("r") = 1.0, py::arg("g_lim") = 0.05, py::arg("step_size") = 1.0, 
        py::arg("max_iter") = 1000);
    
    py::class_<FacialReduction>(m, "FacialReduction")
        .def(py::init<>())
        .def("reduce", &FacialReduction::reduce);
    
    py::class_<res>(m, "res")
        .def_readwrite("sparse_A", &res::sparse_A)
        .def_readwrite("sparse_b", &res::sparse_b)
        .def_readwrite("saved_V", &res::saved_V)
        .def_readwrite("dense_A", &res::dense_A)
        .def_readwrite("dense_b", &res::dense_b)
        .def_readwrite("Q", &res::Q)
        .def_readwrite("z1", &res::z1);
    
    py::class_<SparseRandomWalk, PySparseRandomWalk<>>(m_dense, "SparseRandomWalk")
        .def(py::init<const double >(), py::arg("err") = 0.00001)
        .def("generateCompleteWalk", &SparseRandomWalk::generateCompleteWalk);
    
    py::class_<SparseBallWalk, SparseRandomWalk>(m_sparse, "SparseBallWalk")
        .def(py::init<const double>(), py::arg("r") = 0.3);
    
    py::class_<SparseHitAndRun, SparseRandomWalk>(m_sparse, "SparseHitAndRun")
        .def(py::init<const double, const double>(), 
        py::arg("err") = 0.01, py::arg("r") = 0.5);

    py::class_<SparseBarrierWalk, SparseRandomWalk, PySparseBarrierWalk<>>(m_sparse, "SparseBarrierWalk")
        .def(py::init<const double, const double>(), py::arg("err") = 0.0000001, py::arg("r") = 0.5)
        .def("generateWeight", &SparseBarrierWalk::generateWeight)
        .def("generateCompleteWalk", &SparseRandomWalk::generateCompleteWalk);
    
    py::class_<SparseDikinWalk, SparseBarrierWalk, PySparseBarrierWalk<SparseDikinWalk>>(m_sparse, "SparseDikinWalk")
        .def(py::init<const double, const double>(), py::arg("err") = 0.0000001, py::arg("r") = 0.5);
    
    py::class_<SparseVaidyaWalk, SparseBarrierWalk, PySparseBarrierWalk<SparseVaidyaWalk>>(m_sparse, "SparseVaidyaWalk")
        .def(py::init<const double, const double>(), py::arg("err") = 0.0000001, py::arg("r") = 0.5);
    
    py::class_<SparseJohnWalk, SparseBarrierWalk, PySparseBarrierWalk<SparseJohnWalk>>(m_sparse, "SparseJohnWalk")
        .def(py::init<double, double, double, double, int>(), py::arg("err") = 0.00001, py::arg("r") = 0.5, 
            py::arg("g_lim") = 0.01, py::arg("step_size") = 0.01, py::arg("max_iter") = 1000);
    
    py::class_<SparseDikinLSWalk, SparseBarrierWalk, PySparseBarrierWalk<SparseDikinLSWalk>>(m_sparse, "SparseDikinLSWalk")
        .def(py::init<double, double, double, double, int>(), py::arg("err") = 0.00001, py::arg("r") = 1.2, 
            py::arg("g_lim") = 0.01, py::arg("step_size") = 0.01, py::arg("max_iter") = 1000);
    


}