#ifndef SYSTEM2D_H
#define SYSTEM2D_H

#include <functional>
#include <string>

#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/SparseCore>

#include "utils.h"

using Hamiltonian = Eigen::MatrixXcd;
using SparseHamiltonian = Eigen::SparseMatrix<std::complex<double>>;

struct System2D
{
    std::size_t n_bands = 2;
    std::size_t n_bands_sc = 4;

    // effective mass
    double m = 1.0;

    // lattice constants
    double dx = 1.0;
    double dy = 1.0;

    // Hopping amplitudes
    double tx = 1.0;
    double ty = 1.0;

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

    virtual void set_default() = 0;

    // continues hamiltonians
    virtual Hamiltonian Hk(double kx, double ky) const = 0;
    virtual Hamiltonian HBdG(double kx, double ky) const = 0;

    // continues in kx, discretized in ky hamiltonians
    virtual Hamiltonian HBdG_discrete_ky_onsite(double kx, double y) const = 0;
    virtual Hamiltonian HBdG_discrete_ky_hopping_p(double kx, double y) const = 0;
    virtual Hamiltonian HBdG_discrete_ky_hopping_m(double kx, double y) const = 0;

    virtual Hamiltonian HBdG_discrete_ky(double kx, std::size_t n_ky) const;
    virtual std::vector<Triplet> HBdG_discrete_ky_triplets(double kx, std::size_t n_ky) const;

    // discretized in kx, discretized in ky hamiltonians
    virtual Hamiltonian HBdG_discrete_onsite(double x, double y) const = 0;
    virtual Hamiltonian HBdG_discrete_hopping_xp(double x, double y) const = 0;
    virtual Hamiltonian HBdG_discrete_hopping_xm(double x, double y) const = 0;
    virtual Hamiltonian HBdG_discrete_hopping_yp(double x, double y) const = 0;
    virtual Hamiltonian HBdG_discrete_hopping_ym(double x, double y) const = 0;

    virtual Hamiltonian HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const;
    virtual std::vector<Triplet> HBdG_discrete_triplets(std::size_t n_kx, std::size_t n_ky) const;
};

#endif