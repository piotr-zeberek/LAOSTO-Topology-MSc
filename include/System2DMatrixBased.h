#ifndef SYSTEM2DMATRIXBASED_H
#define SYSTEM2DMATRIXBASED_H

#include <functional>
#include <string>

#include "System2D.h"
#include "utils.h"

struct System2DMatrixBased : public System2D
{
    virtual Hamiltonian HBdG_discrete_ky_onsite(double kx, double y) const = 0;
    virtual Hamiltonian HBdG_discrete_ky_hopping_p(double kx, double y) const = 0;
    virtual Hamiltonian HBdG_discrete_ky_hopping_m(double kx, double y) const = 0;

    virtual Hamiltonian HBdG_discrete_onsite(double x, double y) const = 0;
    virtual Hamiltonian HBdG_discrete_hopping_xp(double x, double y) const = 0;
    virtual Hamiltonian HBdG_discrete_hopping_xm(double x, double y) const = 0;
    virtual Hamiltonian HBdG_discrete_hopping_yp(double x, double y) const = 0;
    virtual Hamiltonian HBdG_discrete_hopping_ym(double x, double y) const = 0;

    virtual Hamiltonian HBdG_discrete_ky(double kx, std::size_t n_ky) const;
    virtual Triplets triplets_HBdG_discrete_ky(double kx, std::size_t n_ky) const;

    virtual Hamiltonian HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const;
    virtual Triplets triplets_HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const;

    Triplets get_triplets(const SparseHamiltonian &H, double tol = 1e-10) const;
    void add_triplets(Triplets &triplets, const Triplets &triplets_to_add, int row_offset, int col_offset) const;

};

#endif