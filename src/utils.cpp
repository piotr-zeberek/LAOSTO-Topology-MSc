#include "utils.h"

#include <pybind11/embed.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace py::literals;

py::scoped_interpreter guard{};
py::module np = py::module::import("numpy");
py::module sp = py::module::import("scipy.sparse");
py::function eigsh = sp.attr("linalg").attr("eigsh");

void adjoint_triplets(std::vector<Triplet> &triplets)
{
    for (auto &t : triplets)
    {
        t = {t.col(), t.row(), std::conj(t.value())};
    }
}

void transpose_triplets(std::vector<Triplet> &triplets)
{
    for (auto &t : triplets)
    {
        t = {t.col(), t.row(), t.value()};
    }
}

void negate_triplets(std::vector<Triplet> &triplets)
{
    for (auto &t : triplets)
    {
        t = {t.row(), t.col(), -t.value()};
    }
}

void negate_transpose_triplets(std::vector<Triplet> &triplets)
{
    for (auto &t : triplets)
    {
        t = {t.col(), t.row(), -t.value()};
    }
}

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

void add_triplets(std::vector<Triplet> &triplets, const std::vector<Triplet> &triplets_to_add, int row_offset, int col_offset)
{
    for (auto &t : triplets_to_add)
    {
        triplets.emplace_back(t.row() + row_offset, t.col() + col_offset, t.value());
    }
}

// // using scipy
// py::object create_sparse_matrix(const std::vector<Triplet> &triplets, std::size_t size)
// {
//     // std::vector<int> rows(triplets.size());
//     // std::vector<int> cols(triplets.size());
//     // std::vector<std::complex<double>> values(triplets.size());

//     // for (int i = 0; i < triplets.size(); ++i)
//     // {
//     //     const auto &t = triplets[i];
//     //     rows[i] = t.row();
//     //     cols[i] = t.col();
//     //     values[i] = t.value();
//     // }

//     // py::array_t<int> rows_array = py::array_t<int>(rows.size(), rows.data());
//     // py::array_t<int> cols_array = py::array_t<int>(cols.size(), cols.data());
//     // py::array_t<std::complex<double>> values_array = py::array_t<std::complex<double>>(values.size(), values.data());

//     py::array_t<int> rows_array = np.attr("zeros")(triplets.size(), np.attr("uint32"));
//     py::array_t<int> cols_array = np.attr("zeros")(triplets.size(), np.attr("uint32"));
//     py::array_t<std::complex<double>> values_array = np.attr("zeros")(triplets.size(), np.attr("complex128"));

//     for (int i = 0; i < triplets.size(); ++i)
//     {
//         const auto &t = triplets[i];
//         rows_array.mutable_at(i) = t.row();
//         cols_array.mutable_at(i) = t.col();
//         values_array.mutable_at(i) = t.value();
//     }

//     py::object sparse_matrix = sp.attr("coo_matrix")(py::make_tuple(values_array, py::make_tuple(rows_array, cols_array)), py::make_tuple(size, size));

//     return sparse_matrix;
// }

Eigen::VectorXd eigenvals_sparse(const std::vector<Triplet> &triplets, std::size_t size, std::size_t n_eigs, double sigma, double tol)
{
    Eigen::SparseMatrix<std::complex<double>> Hsp(size, size);
    Hsp.reserve(triplets.size());
    Hsp.setFromTriplets(triplets.begin(), triplets.end());

    auto kwargs = py::dict("A"_a = Hsp, "k"_a = n_eigs, "sigma"_a = sigma, "ncv"_a = 3 * n_eigs + 10, "tol"_a = tol, "return_eigenvectors"_a = false);

    py::object result = eigsh(**kwargs);

    Eigen::VectorXd eigenvalues = py::cast<Eigen::VectorXd>(result);

    // Sort eigenvalues
    std::sort(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());

    return eigenvalues;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> eigen_sparse(const std::vector<Triplet> &triplets, std::size_t size, std::size_t n_eigs, double sigma, double tol)
{
    Eigen::SparseMatrix<std::complex<double>> Hsp(size, size);
    Hsp.reserve(triplets.size());
    Hsp.setFromTriplets(triplets.begin(), triplets.end());

    auto kwargs = py::dict("A"_a = Hsp, "k"_a = n_eigs, "sigma"_a = sigma, "ncv"_a = 3 * n_eigs + 10, "tol"_a = tol, "return_eigenvectors"_a = true);

    py::object result = eigsh(**kwargs);

    Eigen::VectorXd eigenvalues = py::cast<Eigen::VectorXd>(py::cast<py::tuple>(result)[0]);
    Eigen::MatrixXcd eigenvectors = py::cast<Eigen::MatrixXcd>(py::cast<py::tuple>(result)[1]);

    return {eigenvalues, eigenvectors};
}