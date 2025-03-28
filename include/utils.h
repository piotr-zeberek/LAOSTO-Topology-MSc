#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

// #define ARMA_USE_ARPACK 1
// #define ARMA_USE_SUPERLU 1

#include <armadillo>

using namespace std::complex_literals;
using Triplet = Eigen::Triplet<std::complex<double>>;
using Triplets = std::vector<Triplet>;

void adjoint_triplets(Triplets &triplets);
void transpose_triplets(Triplets &triplets);
void negate_triplets(Triplets &triplets);
void negate_transpose_triplets(Triplets &triplets);

// Pauli matrices
const Eigen::Matrix2cd s0 = Eigen::Matrix2cd::Identity();
const Eigen::Matrix2cd sx{{0, 1}, {1, 0}};
const Eigen::Matrix2cd sy{{0, -1i}, {1i, 0}};
const Eigen::Matrix2cd sz{{1, 0}, {0, -1}};

inline double fmod_positive(double x, double y)
{
    x -= y * static_cast<int>(x / y);
    return x < 0 ? x + y : x;
}

constexpr inline double eV2au(double energy)
{
    return energy * 0.03674932587122423;
}

constexpr inline double meV2au(double energy)
{
    return 0.001 * eV2au(energy);
}

constexpr inline double nm2au(double length)
{
    return length * 18.89726133921252;
}

constexpr inline double T2au(double field)
{
    return field * 4.254382e-6;
}

Eigen::MatrixXcd kron(const Eigen::MatrixXcd &A, const Eigen::MatrixXcd &B);

void add_triplets(Triplets &triplets, const Triplets &triplets_to_add, int row_offset, int col_offset);

Eigen::VectorXd arma_eigenvals_sparse(const Triplets &triplets, std::size_t rows, std::size_t cols, std::size_t n_eigs, double sigma, double tol = meV2au(1e-6));
std::pair<Eigen::VectorXd, Eigen::MatrixXcd> arma_eigen_sparse(const Triplets &triplets, std::size_t rows, std::size_t cols, std::size_t n_eigs, double sigma, double tol = meV2au(1e-6));

#endif