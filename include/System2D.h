#ifndef SYSTEM2D_H
#define SYSTEM2D_H

#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/SparseCore>

#include "utils.h"

using Hamiltonian = Eigen::MatrixXcd;
using HamiltonianFunction = std::function<Hamiltonian(double, double)>;

struct ElementIndex
{
    std::size_t row{};
    std::size_t col{};

    ElementIndex(std::size_t r, std::size_t c) : row(r), col(c) {}
    ElementIndex(const ElementIndex &other) : row(other.row), col(other.col) {}
    ElementIndex() = default;
};

struct System2D
{
    virtual void set_default_parameters() = 0;
    virtual void set_nonzero_indices() = 0;
    virtual void update_HBdG_nonzero_indices();

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

    // full hamiltonians

    // continues hamiltonians
    virtual Hamiltonian Hk(double kx, double ky) const = 0;
    virtual Hamiltonian Delta(double kx, double ky) const = 0;
    virtual Hamiltonian Delta_Adjoint(double kx, double ky) const = 0;
    Hamiltonian mHmkT(double kx, double ky) const
    {
        return -Hk(-kx, -ky).transpose();
    }

    // Hk hamiltonian
    virtual Hamiltonian Hk_discrete_ky(double kx, std::size_t n_ky) const;
    virtual Hamiltonian Hk_discrete(std::size_t n_kx, std::size_t n_ky) const;

    virtual std::vector<Triplet> triplets_Hk_discrete_ky(double kx, std::size_t n_ky) const;
    virtual std::vector<Triplet> triplets_Hk_discrete(std::size_t n_kx, std::size_t n_ky) const;

    // Constructing BdG Hamiltonian
    virtual Hamiltonian HBdG(double kx, double ky) const;
    virtual Hamiltonian HBdG_discrete_ky(double kx, std::size_t n_ky) const;
    virtual Hamiltonian HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const;

    virtual std::vector<Triplet> triplets_HBdG_discrete_ky(double kx, std::size_t n_ky) const;
    virtual std::vector<Triplet> triplets_HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const;

protected:
    // hamiltonian elements

    // discretized in ky
    virtual Hamiltonian Hk_discrete_ky_onsite(double kx, double y) const = 0;
    virtual Hamiltonian Hk_discrete_ky_hopping_p(double kx, double y) const = 0;
    virtual Hamiltonian Hk_discrete_ky_hopping_m(double kx, double y) const = 0;

    virtual Hamiltonian Delta_discrete_ky_onsite(double kx, double y) const = 0;
    virtual Hamiltonian Delta_discrete_ky_hopping_p(double kx, double y) const = 0;
    virtual Hamiltonian Delta_discrete_ky_hopping_m(double kx, double y) const = 0;

    virtual Hamiltonian Delta_Adjoint_discrete_ky_onsite(double kx, double y) const = 0;
    virtual Hamiltonian Delta_Adjoint_discrete_ky_hopping_p(double kx, double y) const = 0;
    virtual Hamiltonian Delta_Adjoint_discrete_ky_hopping_m(double kx, double y) const = 0;

    virtual Hamiltonian mHmkT_discrete_ky_onsite(double kx, double y) const = 0;
    virtual Hamiltonian mHmkT_discrete_ky_hopping_p(double kx, double y) const = 0;
    virtual Hamiltonian mHmkT_discrete_ky_hopping_m(double kx, double y) const = 0;

    virtual Hamiltonian HBdG_discrete_ky_onsite(double kx, double y) const;
    virtual Hamiltonian HBdG_discrete_ky_hopping_p(double kx, double y) const;
    virtual Hamiltonian HBdG_discrete_ky_hopping_m(double kx, double y) const;

    // discretized in kx,ky
    virtual Hamiltonian Hk_discrete_onsite(double x, double y) const = 0;
    virtual Hamiltonian Hk_discrete_hopping_xp(double x, double y) const = 0;
    virtual Hamiltonian Hk_discrete_hopping_xm(double x, double y) const = 0;
    virtual Hamiltonian Hk_discrete_hopping_yp(double x, double y) const = 0;
    virtual Hamiltonian Hk_discrete_hopping_ym(double x, double y) const = 0;

    virtual Hamiltonian Delta_discrete_onsite(double x, double y) const = 0;
    virtual Hamiltonian Delta_discrete_hopping_xp(double x, double y) const = 0;
    virtual Hamiltonian Delta_discrete_hopping_xm(double x, double y) const = 0;
    virtual Hamiltonian Delta_discrete_hopping_yp(double x, double y) const = 0;
    virtual Hamiltonian Delta_discrete_hopping_ym(double x, double y) const = 0;

    virtual Hamiltonian Delta_Adjoint_discrete_onsite(double x, double y) const = 0;
    virtual Hamiltonian Delta_Adjoint_discrete_hopping_xp(double x, double y) const = 0;
    virtual Hamiltonian Delta_Adjoint_discrete_hopping_xm(double x, double y) const = 0;
    virtual Hamiltonian Delta_Adjoint_discrete_hopping_yp(double x, double y) const = 0;
    virtual Hamiltonian Delta_Adjoint_discrete_hopping_ym(double x, double y) const = 0;

    virtual Hamiltonian mHmkT_discrete_onsite(double x, double y) const = 0;
    virtual Hamiltonian mHmkT_discrete_hopping_xp(double x, double y) const = 0;
    virtual Hamiltonian mHmkT_discrete_hopping_xm(double x, double y) const = 0;
    virtual Hamiltonian mHmkT_discrete_hopping_yp(double x, double y) const = 0;
    virtual Hamiltonian mHmkT_discrete_hopping_ym(double x, double y) const = 0;

    virtual Hamiltonian HBdG_discrete_onsite(double x, double y) const;
    virtual Hamiltonian HBdG_discrete_hopping_xp(double x, double y) const;
    virtual Hamiltonian HBdG_discrete_hopping_xm(double x, double y) const;
    virtual Hamiltonian HBdG_discrete_hopping_yp(double x, double y) const;
    virtual Hamiltonian HBdG_discrete_hopping_ym(double x, double y) const;

    // nonzero indices for triplets construction
    std::vector<ElementIndex> Hk_discrete_ky_onsite_nonzero_indices;
    std::vector<ElementIndex> Hk_discrete_ky_hopping_p_nonzero_indices;
    std::vector<ElementIndex> Hk_discrete_ky_hopping_m_nonzero_indices;

    std::vector<ElementIndex> Delta_discrete_ky_onsite_nonzero_indices;
    std::vector<ElementIndex> Delta_discrete_ky_hopping_p_nonzero_indices;
    std::vector<ElementIndex> Delta_discrete_ky_hopping_m_nonzero_indices;

    std::vector<ElementIndex> Delta_Adjoint_discrete_ky_onsite_nonzero_indices;
    std::vector<ElementIndex> Delta_Adjoint_discrete_ky_hopping_p_nonzero_indices;
    std::vector<ElementIndex> Delta_Adjoint_discrete_ky_hopping_m_nonzero_indices;

    std::vector<ElementIndex> mHmkT_discrete_ky_onsite_nonzero_indices;
    std::vector<ElementIndex> mHmkT_discrete_ky_hopping_p_nonzero_indices;
    std::vector<ElementIndex> mHmkT_discrete_ky_hopping_m_nonzero_indices;

    std::vector<ElementIndex> HBdG_discrete_ky_onsite_nonzero_indices;
    std::vector<ElementIndex> HBdG_discrete_ky_hopping_p_nonzero_indices;
    std::vector<ElementIndex> HBdG_discrete_ky_hopping_m_nonzero_indices;

    std::vector<ElementIndex> Hk_discrete_onsite_nonzero_indices;
    std::vector<ElementIndex> Hk_discrete_hopping_xp_nonzero_indices;
    std::vector<ElementIndex> Hk_discrete_hopping_xm_nonzero_indices;
    std::vector<ElementIndex> Hk_discrete_hopping_yp_nonzero_indices;
    std::vector<ElementIndex> Hk_discrete_hopping_ym_nonzero_indices;

    std::vector<ElementIndex> Delta_discrete_onsite_nonzero_indices;
    std::vector<ElementIndex> Delta_discrete_hopping_xp_nonzero_indices;
    std::vector<ElementIndex> Delta_discrete_hopping_xm_nonzero_indices;
    std::vector<ElementIndex> Delta_discrete_hopping_yp_nonzero_indices;
    std::vector<ElementIndex> Delta_discrete_hopping_ym_nonzero_indices;

    std::vector<ElementIndex> Delta_Adjoint_discrete_onsite_nonzero_indices;
    std::vector<ElementIndex> Delta_Adjoint_discrete_hopping_xp_nonzero_indices;
    std::vector<ElementIndex> Delta_Adjoint_discrete_hopping_xm_nonzero_indices;
    std::vector<ElementIndex> Delta_Adjoint_discrete_hopping_yp_nonzero_indices;
    std::vector<ElementIndex> Delta_Adjoint_discrete_hopping_ym_nonzero_indices;

    std::vector<ElementIndex> mHmkT_discrete_onsite_nonzero_indices;
    std::vector<ElementIndex> mHmkT_discrete_hopping_xp_nonzero_indices;
    std::vector<ElementIndex> mHmkT_discrete_hopping_xm_nonzero_indices;
    std::vector<ElementIndex> mHmkT_discrete_hopping_yp_nonzero_indices;
    std::vector<ElementIndex> mHmkT_discrete_hopping_ym_nonzero_indices;

    std::vector<ElementIndex> HBdG_discrete_onsite_nonzero_indices;
    std::vector<ElementIndex> HBdG_discrete_hopping_xp_nonzero_indices;
    std::vector<ElementIndex> HBdG_discrete_hopping_xm_nonzero_indices;
    std::vector<ElementIndex> HBdG_discrete_hopping_yp_nonzero_indices;
    std::vector<ElementIndex> HBdG_discrete_hopping_ym_nonzero_indices;

    Hamiltonian assemble_HBdG(const Hamiltonian &Hk_mat, const Hamiltonian &Delta_mat, const Hamiltonian &Delta_Adjoint_mat, const Hamiltonian &mHmkT_mat) const;
    std::vector<ElementIndex> assemble_HBdG_nonzero_indices(const std::vector<ElementIndex> &Hk_nonzero_indices,
                                                        const std::vector<ElementIndex> &Delta_nonzero_indices,
                                                        const std::vector<ElementIndex> &Delta_Adjoint_nonzero_indices,
                                                        const std::vector<ElementIndex> &mHmkT_nonzero_indices) const;
    std::vector<Triplet> assemble_triplets_HBdG(const std::vector<Triplet> &Hk_triplets, const std::vector<Triplet> &Delta_triplets, const std::vector<Triplet> &Delta_adjoint_triplets, const std::vector<Triplet> &mHmkT_triplets) const;

    Hamiltonian assemble_matrix_discrete_ky(double kx, std::size_t n_ky,
                                            const HamiltonianFunction &onsite,
                                            const HamiltonianFunction &hopping_p,
                                            const HamiltonianFunction &hopping_m) const;

    Hamiltonian assemble_matrix_discrete(std::size_t n_kx, std::size_t n_ky,
                                         const HamiltonianFunction &onsite,
                                         const HamiltonianFunction &hopping_xp,
                                         const HamiltonianFunction &hopping_xm,
                                         const HamiltonianFunction &hopping_yp,
                                         const HamiltonianFunction &hopping_ym) const;

    std::vector<Triplet> assemble_triplets_discrete_ky(double kx, std::size_t n_ky,
                                                       const HamiltonianFunction &onsite,
                                                       const HamiltonianFunction &hopping_p,
                                                       const HamiltonianFunction &hopping_m,
                                                       const std::vector<ElementIndex> &onsite_nonzero_indices,
                                                       const std::vector<ElementIndex> &hopping_p_nonzero_indices,
                                                       const std::vector<ElementIndex> &hopping_m_nonzero_indices) const;

    std::vector<Triplet> assemble_triplets_discrete(std::size_t n_kx, std::size_t n_ky,
                                                    const HamiltonianFunction &onsite,
                                                    const HamiltonianFunction &hopping_xp,
                                                    const HamiltonianFunction &hopping_xm,
                                                    const HamiltonianFunction &hopping_yp,
                                                    const HamiltonianFunction &hopping_ym,
                                                    const std::vector<ElementIndex> &onsite_nonzero_indices,
                                                    const std::vector<ElementIndex> &hopping_xp_nonzero_indices,
                                                    const std::vector<ElementIndex> &hopping_xm_nonzero_indices,
                                                    const std::vector<ElementIndex> &hopping_yp_nonzero_indices,
                                                    const std::vector<ElementIndex> &hopping_ym_nonzero_indices) const;
};

#endif // SYSTEM2D_H