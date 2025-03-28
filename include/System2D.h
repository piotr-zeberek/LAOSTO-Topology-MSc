#ifndef SYSTEM2D_H
#define SYSTEM2D_H

#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/SparseCore>

#include "utils.h"

using Hamiltonian = Eigen::MatrixXcd;
using SparseHamiltonian = Eigen::SparseMatrix<std::complex<double>>;

struct System2D
{
    virtual void set_default_parameters() = 0;

    std::size_t n_bands = 2;

    // lattice constants
    double dx = 1.0;
    double dy = 1.0;

    // Chemical potential
    double mu{};

    // Lande g-factor
    double g_Lande{};

    // External magnetic field
    double Bx{};
    double By{};
    double Bz{};

    // Superconducting energy gap
    double delta_SC{};

    // continues hamiltonians
    virtual Hamiltonian Hk(double kx, double ky) const = 0;
    virtual Hamiltonian HBdG(double kx, double ky) const = 0;
    virtual Hamiltonian HBdG_discrete_ky(double kx, std::size_t n_ky) const = 0;
    virtual Hamiltonian HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const = 0;
    
    virtual std::vector<Triplet> triplets_HBdG_discrete_ky(double kx, std::size_t n_ky) const = 0;
    virtual std::vector<Triplet> triplets_HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const = 0;
};

#endif // SYSTEM2D_H