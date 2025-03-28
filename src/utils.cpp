#include "utils.h"

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

Eigen::VectorXd arma_eigenvals_sparse(const std::vector<Triplet> &triplets, std::size_t rows, std::size_t cols, std::size_t n_eigs, double sigma, double tol)
{
    arma::sp_cx_mat armaHsp(rows, cols);

    for (const auto &t : triplets)
    {
        armaHsp(t.row(), t.col()) = t.value();
    }

    arma::eigs_opts opts;
    opts.tol = tol;
    opts.subdim = 3 * n_eigs + 10;

    arma::cx_vec eigval;
    arma::eigs_gen(eigval, armaHsp, n_eigs, sigma, opts);

    Eigen::VectorXd evals(n_eigs);

    for (auto i = 0; i < eigval.size(); ++i)
    {
        evals(i) = eigval(i).real();
    }

    return evals;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> arma_eigen_sparse(const std::vector<Triplet> &triplets, std::size_t rows, std::size_t cols, std::size_t n_eigs, double sigma, double tol)
{
    arma::sp_cx_mat armaHsp(rows, cols);

    for (const auto &t : triplets)
    {
        armaHsp(t.row(), t.col()) = t.value();
    }

    arma::eigs_opts opts;
    opts.tol = tol;
    opts.subdim = 3 * n_eigs + 10;

    arma::cx_vec eigval;
    arma::cx_mat eigvec;
    arma::eigs_gen(eigval, eigvec, armaHsp, n_eigs, sigma, opts);

    Eigen::VectorXd evals(n_eigs);
    Eigen::MatrixXcd evecs(rows, n_eigs);

    for (auto i = 0; i < eigval.size(); ++i)
    {
        evals(i) = eigval(i).real();

        for (auto j = 0; j < rows; ++j)
        {
            evecs(j, i) = eigvec(j, i);
        }
    }

    return {evals, evecs};
}