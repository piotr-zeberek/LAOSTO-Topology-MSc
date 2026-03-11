#include "utils.h"
#include <iostream>

#include <pybind11/embed.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace py::literals;

py::scoped_interpreter guard{};
py::module np = py::module::import("numpy");
py::module sp = py::module::import("scipy.sparse");
py::function eigsh = sp.attr("linalg").attr("eigsh");

py::module pf = py::module::import("pfapack.ctypes");
py::function cpf = pf.attr("pfaffian");

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

Eigen::VectorXd orbital_prob_den(const Eigen::VectorXcd &eigenvec, std::size_t n_orbitals)
{
    Eigen::VectorXd prob_den = Eigen::VectorXd::Zero(n_orbitals);

    std::size_t n_sites = eigenvec.size() / n_orbitals;

    for (std::size_t i = 0; i < n_sites; ++i)
    {
        prob_den += eigenvec.segment(i * n_orbitals, n_orbitals).cwiseAbs2();
    }

    return prob_den;
}

Eigen::VectorXd eigenvals_sparse(const SparseHamiltonian &sparse_H, std::size_t n_eigs, double sigma, double tol)
{
    auto kwargs = py::dict("A"_a = sparse_H, "k"_a = n_eigs, "sigma"_a = sigma, "ncv"_a = 3 * n_eigs + 10, "tol"_a = tol, "return_eigenvectors"_a = false);

    py::object result = eigsh(**kwargs);

    Eigen::VectorXd eigenvalues = py::cast<Eigen::VectorXd>(result);

    // Sort eigenvalues
    std::sort(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());

    return eigenvalues;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> eigen_sparse(const SparseHamiltonian &sparse_H, std::size_t n_eigs, double sigma, double tol)
{
    auto kwargs = py::dict("A"_a = sparse_H, "k"_a = n_eigs, "sigma"_a = sigma, "ncv"_a = 3 * n_eigs + 10, "tol"_a = tol, "return_eigenvectors"_a = true);

    py::object result = eigsh(**kwargs);

    Eigen::VectorXd eigenvalues = py::cast<Eigen::VectorXd>(py::cast<py::tuple>(result)[0]);
    Eigen::MatrixXcd eigenvectors = py::cast<Eigen::MatrixXcd>(py::cast<py::tuple>(result)[1]);

    // Sort eigenvalues and eigenvectors
    std::vector<std::size_t> indices(eigenvalues.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&eigenvalues](std::size_t i1, std::size_t i2)
              { return eigenvalues(i1) < eigenvalues(i2); });

    Eigen::VectorXd sorted_eigenvalues(eigenvalues.size());
    Eigen::MatrixXcd sorted_eigenvectors(eigenvectors.rows(), eigenvectors.cols());
    for (std::size_t i = 0; i < indices.size(); ++i)
    {
        sorted_eigenvalues(i) = eigenvalues(indices[i]);
        sorted_eigenvectors.col(i) = eigenvectors.col(indices[i]);
    }

    return {sorted_eigenvalues, sorted_eigenvectors};
}

double pfaffian(const Eigen::MatrixXd &A)
{
    // using pfapack C extension
    auto kwargs = py::dict("matrix"_a = A / A.cwiseAbs().maxCoeff(), "method"_a = "P", "avoid_overflow"_a = true);

    return cpf(**kwargs).cast<double>();
}