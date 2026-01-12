#ifndef LAOSTO_H
#define LAOSTO_H

#include "System2D.h"
#include "utils.h"

#include <unsupported/Eigen/KroneckerProduct>

struct LAOSTO : public System2D
{
    // Hopping parameters: l - light, h - heavy, d - coupling between the dxz /dyz
    double tl{};
    double th{};
    double td{};

    // Energy difference between xy and xz/yz orbitals
    double delta_E{};

    // Atomic spin-orbit coupling
    double delta_SO{};

    // Rashba spin-orbit coupling
    double delta_RSO{};

    LAOSTO() : System2D()
    {
        set_default_parameters();

        // Atomic SO matrix
        Hamiltonian H(n_bands, n_bands);
        Eigen::Matrix2cd z = Eigen::Matrix2cd::Zero();
        H << z, 1i * sx, -1i * sy,
            -1i * sx, z, 1i * sz,
            1i * sy, -1i * sz, z;

        HSO_mat = H;
    }

    void set_default_parameters() override
    {
        n_bands = 6;

        double dx_nm = 0.39;
        double dy_nm = 0.39;

        dx = nm2au(dx_nm);
        dy = nm2au(dy_nm);

        tl = meV2au(875.0);
        th = meV2au(40.0);
        td = meV2au(40.0);

        delta_E = meV2au(47.0);
        delta_SO = meV2au(10.0);
        delta_RSO = meV2au(20.0);
        delta_SC = meV2au(0.02);

        g_Lande = 3.0;
        Bx = 0.0;
        By = 0.0;
        Bz = 0.0;
    }

protected:
    // Triplets only upper triangular part of the entire hamiltonian for onsites
    // continues hamiltonians
    std::vector<Triplet> Hk_triplets(double kx, double ky) const;
    std::vector<Triplet> Delta_triplets(double kx, double ky) const;
    std::vector<Triplet> mHmkT_triplets(double kx, double ky) const;

    // discrete in y direction
    std::vector<Triplet> Hk_discrete_ky_onsite_triplets(double kx, double y) const;
    std::vector<Triplet> Hk_discrete_ky_hopping_p_triplets(double kx, double y) const;

    std::vector<Triplet> Delta_discrete_ky_onsite_triplets(double kx, double y) const;
    // std::vector<Triplet> Delta_discrete_ky_hopping_p_triplets(double kx, double y) const;

    std::vector<Triplet> mHmkT_discrete_ky_onsite_triplets(double kx, double y) const;
    std::vector<Triplet> mHmkT_discrete_ky_hopping_p_triplets(double kx, double y) const;

    // discrete in both x and y directions
    std::vector<Triplet> Hk_discrete_onsite_triplets(double x, double y) const;
    std::vector<Triplet> Hk_discrete_hopping_xp_triplets(double x, double y) const;
    std::vector<Triplet> Hk_discrete_hopping_yp_triplets(double x, double y) const;
    std::vector<Triplet> Hk_discrete_hopping_pp_triplets(double x, double y) const;
    std::vector<Triplet> Hk_discrete_hopping_pm_triplets(double x, double y) const;

    std::vector<Triplet> Delta_discrete_onsite_triplets(double x, double y) const;
    // std::vector<Triplet> Delta_discrete_hopping_xp_triplets(double x, double y) const;
    // std::vector<Triplet> Delta_discrete_hopping_yp_triplets(double x, double y) const;
    // std::vector<Triplet> Delta_discrete_hopping_pp_triplets(double x, double y) const;
    // std::vector<Triplet> Delta_discrete_hopping_pm_triplets(double x, double y) const;

    std::vector<Triplet> mHmkT_discrete_onsite_triplets(double x, double y) const;
    std::vector<Triplet> mHmkT_discrete_hopping_xp_triplets(double x, double y) const;
    std::vector<Triplet> mHmkT_discrete_hopping_yp_triplets(double x, double y) const;
    std::vector<Triplet> mHmkT_discrete_hopping_pp_triplets(double x, double y) const;
    std::vector<Triplet> mHmkT_discrete_hopping_pm_triplets(double x, double y) const;

public:
    // Atomic SO matrix
    Hamiltonian HSO_mat;

    // angular momentum matrices
    const Eigen::Matrix3cd Lx{{0, 1i, 0},
                              {-1i, 0, 0},
                              {0, 0, 0}};
    const Eigen::Matrix3cd Ly{{0, 0, -1i},
                              {0, 0, 0},
                              {1i, 0, 0}};
    const Eigen::Matrix3cd Lz{{0, 0, 0},
                              {0, 0, 1i},
                              {0, -1i, 0}};

    // 3x3 complex identity matrix
    const Eigen::Matrix3cd I3 = Eigen::Matrix3cd::Identity();

    // for HZeeman
    const Hamiltonian Lxs0 = Eigen::kroneckerProduct(Lx, s0);
    const Hamiltonian Lys0 = Eigen::kroneckerProduct(Ly, s0);
    const Hamiltonian Lzs0 = Eigen::kroneckerProduct(Lz, s0);
    const Hamiltonian I3sx = Eigen::kroneckerProduct(I3, sx);
    const Hamiltonian I3sy = Eigen::kroneckerProduct(I3, sy);
    const Hamiltonian I3sz = Eigen::kroneckerProduct(I3, sz);

    Hamiltonian Hk_mat(double kx, double ky) const;
    Hamiltonian Delta(double kx, double ky) const;

    // discretized in ky
    Hamiltonian Hk_discrete_ky_onsite(double kx, double y) const;
    Hamiltonian Hk_discrete_ky_hopping_p(double kx, double y) const;

    Hamiltonian mHmkT_discrete_ky_onsite(double kx, double y) const;
    Hamiltonian mHmkT_discrete_ky_hopping_p(double kx, double y) const;

    // discretized in kx,ky
    Hamiltonian Hk_discrete_onsite(double x, double y) const;
    Hamiltonian Hk_discrete_hopping_xp(double x, double y) const;
    Hamiltonian Hk_discrete_hopping_yp(double x, double y) const;

    Hamiltonian mHmkT_discrete_onsite(double x, double y) const;
    Hamiltonian mHmkT_discrete_hopping_xp(double x, double y) const;
    Hamiltonian mHmkT_discrete_hopping_yp(double x, double y) const;

    // additional diagonal hoppings
    Hamiltonian Hk_discrete_hopping_pp(double x, double y) const;
    Hamiltonian Hk_discrete_hopping_pm(double x, double y) const;

    Hamiltonian mHmkT_discrete_hopping_pp(double x, double y) const;
    Hamiltonian mHmkT_discrete_hopping_pm(double x, double y) const;

    Hamiltonian Hkin(double kx, double ky) const;
    Hamiltonian Hkin(double kx) const;
    Hamiltonian Hkin() const;
    Hamiltonian HZeeman() const;
    Hamiltonian HAtomicSO() const;
    Hamiltonian HRashba(double kx, double ky) const;
    Hamiltonian HRashba(double kx) const;

    double Ekxy(double kx, double ky) const
    {
        return 2.0 * tl * (2.0 - std::cos(kx) - std::cos(ky)) - delta_E;
    }

    double Ekxz(double kx, double ky) const
    {
        return 2.0 * tl * (1.0 - std::cos(kx)) + 2.0 * th * (1.0 - std::cos(ky));
    }
    double Ekyz(double kx, double ky) const
    {
        return 2.0 * tl * (1.0 - std::cos(ky)) + 2.0 * th * (1.0 - std::cos(kx));
    }

    double Ek_h(double kx, double ky) const
    {
        return 2.0 * td * std::sin(kx) * std::sin(ky);
    }

    double Ekxy(double kx) const
    {
        return 2.0 * tl * (2.0 - std::cos(kx)) - delta_E;
    }

    double Ekxz(double kx) const
    {
        return 2.0 * tl * (1.0 - std::cos(kx)) + 2.0 * th;
    }
    double Ekyz(double kx) const
    {
        return 2.0 * tl + 2.0 * th * (1.0 - std::cos(kx));
    }

    std::complex<double> Ek_h(double kx) const
    {
        return -1i * td * std::sin(kx);
    }

    double Ekxy() const
    {
        return 4.0 * tl - delta_E;
    }

    double Ekxz() const
    {
        return 2.0 * (tl + th);
    }

    double Ekyz() const
    {
        return 2.0 * (tl + th);
    }
};

#endif