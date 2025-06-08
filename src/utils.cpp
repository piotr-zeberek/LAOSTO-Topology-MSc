#include "utils.h"

#include <pybind11/embed.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace py::literals;

py::scoped_interpreter guard{};
py::module np = py::module::import("numpy");
py::module sp = py::module::import("scipy.sparse");
py::function eigsh = sp.attr("linalg").attr("eigsh");


Eigen::MatrixXcd kron(const Eigen::MatrixXcd &A, const Eigen::MatrixXcd &B)
{
    Eigen::MatrixXcd res(A.rows() * B.rows(), A.cols() * B.cols());

    for (auto i = 0; i < A.rows(); ++i)
    {
        for (auto j = 0; j < A.cols(); ++j)
        {
            res.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i, j) * B;
        }
    }

    return res;
}


Eigen::VectorXd eigenvals_sparse(const SparseHamiltonian& sparse_H, std::size_t n_eigs, double sigma, double tol)
{
    auto kwargs = py::dict("A"_a = sparse_H, "k"_a = n_eigs, "sigma"_a = sigma, "ncv"_a = 3 * n_eigs + 10, "tol"_a = tol, "return_eigenvectors"_a = false);

    py::object result = eigsh(**kwargs);

    Eigen::VectorXd eigenvalues = py::cast<Eigen::VectorXd>(result);

    // Sort eigenvalues
    std::sort(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());

    return eigenvalues;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> eigen_sparse(const SparseHamiltonian& sparse_H, std::size_t n_eigs, double sigma, double tol)
{
    auto kwargs = py::dict("A"_a = sparse_H, "k"_a = n_eigs, "sigma"_a = sigma, "ncv"_a = 3 * n_eigs + 10, "tol"_a = tol, "return_eigenvectors"_a = true);

    py::object result = eigsh(**kwargs);

    Eigen::VectorXd eigenvalues = py::cast<Eigen::VectorXd>(py::cast<py::tuple>(result)[0]);
    Eigen::MatrixXcd eigenvectors = py::cast<Eigen::MatrixXcd>(py::cast<py::tuple>(result)[1]);

    return {eigenvalues, eigenvectors};
}