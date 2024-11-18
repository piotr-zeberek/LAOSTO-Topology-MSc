#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Core>

using namespace std::complex_literals;

// Pauli matrices
const Eigen::Matrix2cd s0 = Eigen::Matrix2cd::Identity();
const Eigen::Matrix2cd sx{{0, 1}, {1, 0}};
const Eigen::Matrix2cd sy{{0, -1i}, {1i, 0}};
const Eigen::Matrix2cd sz{{1, 0}, {0, -1}};

// angular momentum matrices
const Eigen::Matrix3cd Lx{{0, 1i, 0}, {-1i, 0, 0}, {0, 0, 0}};
const Eigen::Matrix3cd Ly{{0, 0, -1i}, {0, 0, 0}, {1i, 0, 0}};
const Eigen::Matrix3cd Lz{{0, 0, 0}, {0, -1i, 0}, {0, 0, 1i}};

Eigen::MatrixXcd kron(const Eigen::MatrixXcd &A, const Eigen::MatrixXcd &B);

inline double fmod_positive(double x, double y){
    x-=y*static_cast<int>(x/y);
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

#endif