#ifndef SYSTEM2D_H
#define SYSTEM2D_H

#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/SparseCore>

#include "utils.h"

struct System2D
{
    // in normal state
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

    virtual ~System2D() = default;

    virtual void set_default_parameters() = 0;

    // Normal Hamiltonian
    Hamiltonian Hk(double kx, double ky) const;

    Hamiltonian Hk_discrete_ky(double kx, std::size_t n_ky) const;
    Hamiltonian Hk_discrete(std::size_t n_kx, std::size_t n_ky) const;

    SparseHamiltonian Hk_discrete_ky_sparse(double kx, std::size_t n_ky) const;
    SparseHamiltonian Hk_discrete_sparse(std::size_t n_kx, std::size_t n_ky) const;

    // BdG Hamiltonian
    Hamiltonian HBdG(double kx, double ky) const;

    Hamiltonian HBdG_discrete_ky(double kx, std::size_t n_ky) const;
    Hamiltonian HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const;

    SparseHamiltonian HBdG_discrete_ky_sparse(double kx, std::size_t n_ky) const;
    SparseHamiltonian HBdG_discrete_sparse(std::size_t n_kx, std::size_t n_ky) const;

protected:

    std::vector<Triplet> generate_triplets(const Hamiltonian &H, bool only_upper_triangular = false) const;
    std::vector<Triplet> negate_triplets(const std::vector<Triplet> &triplets) const;
    // Triplets only upper triangular part of the entire hamiltonian for onsites
    // continues hamiltonians
    virtual std::vector<Triplet> Hk_triplets(double kx, double ky) const { return {}; }
    virtual std::vector<Triplet> Delta_triplets(double kx, double ky) const { return {}; }
    virtual std::vector<Triplet> mHmkT_triplets(double kx, double ky) const { return {}; }

    // discrete in y direction
    virtual std::vector<Triplet> Hk_discrete_ky_onsite_triplets(double kx, double y) const { return {}; }
    virtual std::vector<Triplet> Hk_discrete_ky_hopping_p_triplets(double kx, double y) const { return {}; }

    virtual std::vector<Triplet> Delta_discrete_ky_onsite_triplets(double kx, double y) const { return {}; }
    virtual std::vector<Triplet> Delta_discrete_ky_hopping_p_triplets(double kx, double y) const { return {}; }

    virtual std::vector<Triplet> mHmkT_discrete_ky_onsite_triplets(double kx, double y) const { return {}; }
    virtual std::vector<Triplet> mHmkT_discrete_ky_hopping_p_triplets(double kx, double y) const { return {}; }

    // discrete in both x and y directions
    virtual std::vector<Triplet> Hk_discrete_onsite_triplets(double x, double y) const { return {}; }
    virtual std::vector<Triplet> Hk_discrete_hopping_xp_triplets(double x, double y) const { return {}; }
    virtual std::vector<Triplet> Hk_discrete_hopping_yp_triplets(double x, double y) const { return {}; }
    virtual std::vector<Triplet> Hk_discrete_hopping_pp_triplets(double x, double y) const { return {}; }
    virtual std::vector<Triplet> Hk_discrete_hopping_pm_triplets(double x, double y) const { return {}; }

    virtual std::vector<Triplet> Delta_discrete_onsite_triplets(double x, double y) const { return {}; }
    virtual std::vector<Triplet> Delta_discrete_hopping_xp_triplets(double x, double y) const { return {}; }
    virtual std::vector<Triplet> Delta_discrete_hopping_yp_triplets(double x, double y) const { return {}; }
    virtual std::vector<Triplet> Delta_discrete_hopping_pp_triplets(double x, double y) const { return {}; }
    virtual std::vector<Triplet> Delta_discrete_hopping_pm_triplets(double x, double y) const { return {}; }

    virtual std::vector<Triplet> mHmkT_discrete_onsite_triplets(double x, double y) const { return {}; }
    virtual std::vector<Triplet> mHmkT_discrete_hopping_xp_triplets(double x, double y) const { return {}; }
    virtual std::vector<Triplet> mHmkT_discrete_hopping_yp_triplets(double x, double y) const { return {}; }
    virtual std::vector<Triplet> mHmkT_discrete_hopping_pp_triplets(double x, double y) const { return {}; }
    virtual std::vector<Triplet> mHmkT_discrete_hopping_pm_triplets(double x, double y) const { return {}; }

    // // Triplets for the normal Hamiltonian
private:
    std::vector<Triplet> Hk_discrete_ky_triplets(double kx, std::size_t n_ky) const;
    std::vector<Triplet> Hk_discrete_triplets(std::size_t n_kx, std::size_t n_ky) const;

    // // Triplets for the BdG Hamiltonian
    std::vector<Triplet> join_triplets_for_HBdG(const std::vector<Triplet> &Hk_tr,
                                                const std::vector<Triplet> &Delta_tr,
                                                const std::vector<Triplet> &mHmkT_tr) const;
    std::vector<Triplet> HBdG_triplets(double kx, double ky) const;
    std::vector<Triplet> HBdG_discrete_ky_triplets(double kx, std::size_t n_ky) const;
    std::vector<Triplet> HBdG_discrete_triplets(std::size_t n_kx, std::size_t n_ky) const;

    void append_triplets(std::size_t row_offset,
                         std::size_t col_offset,
                         std::vector<Triplet> &target,
                         const std::vector<Triplet> &source) const;

    std::vector<Triplet> assemble_triplets_discrete_ky(const TripletFunc &onsite_tf,
                                                       const TripletFunc &hopping_p_tf,
                                                       double kx, std::size_t n_ky, std::size_t submatrix_size) const;

    std::vector<Triplet> assemble_triplets_discrete(const TripletFunc &onsite_tf,
                                                    const TripletFunc &hopping_xp_tf,
                                                    const TripletFunc &hopping_yp_tf,
                                                    const TripletFunc &hopping_pp_tf,
                                                    const TripletFunc &hopping_pm_tf,
                                                    std::size_t n_kx, std::size_t n_ky, std::size_t submatrix_size) const;

    Hamiltonian assemble_matrix(const std::vector<Triplet> &triplets, std::size_t n_rows, std::size_t n_cols) const;
    SparseHamiltonian assemble_sparse_matrix(const std::vector<Triplet> &triplets, std::size_t n_rows, std::size_t n_cols) const;
};

#endif // SYSTEM2D_H