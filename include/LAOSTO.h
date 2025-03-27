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

    // angular momentum matrices
    const Eigen::Matrix3cd Lx{{0, 1i, 0}, {-1i, 0, 0}, {0, 0, 0}};
    const Eigen::Matrix3cd Ly{{0, 0, -1i}, {0, 0, 0}, {1i, 0, 0}};
    const Eigen::Matrix3cd Lz{{0, 0, 0}, {0, -1i, 0}, {0, 0, 1i}};

    // 3x3 complex identity matrix
    const Eigen::Matrix3cd I3 = Eigen::Matrix3cd::Identity();

    // for HZeeman
    const Hamiltonian Lxs0 = Eigen::kroneckerProduct(Lx, s0);
    const Hamiltonian Lys0 = Eigen::kroneckerProduct(Ly, s0);
    const Hamiltonian Lzs0 = Eigen::kroneckerProduct(Lz, s0);
    const Hamiltonian I3sx = Eigen::kroneckerProduct(I3, sx);
    const Hamiltonian I3sy = Eigen::kroneckerProduct(I3, sy);
    const Hamiltonian I3sz = Eigen::kroneckerProduct(I3, sz);

    // for SC
    const Hamiltonian syI3sy = Eigen::kroneckerProduct(sy, Eigen::kroneckerProduct(I3, sy));

    // Atomic SO matrix
    Hamiltonian HSO_mat;

    LAOSTO()
    {
        set_default();

        // Atomic SO matrix
        Hamiltonian H(n_bands, n_bands);
        Eigen::Matrix2cd z = Eigen::Matrix2cd::Zero();
        H << z, 1i * sx, -1i * sy,
            -1i * sx, z, 1i * sz,
            1i * sy, -1i * sz, z;

        HSO_mat = H;
    }

    void set_default() override
    {
        n_bands = 6;
        n_bands_sc = 2 * n_bands;

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
        delta_SC = meV2au(0.2);

        mu = 0.0;

        g_Lande = 3.0;
        Bx = 0.0;
        By = 0.0;
        Bz = 0.0;
    }

    Hamiltonian Hk(double kx, double ky) const override;
    Hamiltonian HBdG(double kx, double ky) const override;

    Hamiltonian HBdG_discrete_ky_onsite(double kx, double y) const override;
    Hamiltonian HBdG_discrete_ky_hopping_p(double kx, double y) const override;
    Hamiltonian HBdG_discrete_ky_hopping_m(double kx, double y) const override;

    Hamiltonian HBdG_discrete_onsite(double x, double y) const override;
    Hamiltonian HBdG_discrete_hopping_xp(double x, double y) const override;
    Hamiltonian HBdG_discrete_hopping_xm(double x, double y) const override;
    Hamiltonian HBdG_discrete_hopping_yp(double x, double y) const override;
    Hamiltonian HBdG_discrete_hopping_ym(double x, double y) const override;

    Hamiltonian Hkin(double kx, double ky) const;
    Hamiltonian HZeeman() const;
    Hamiltonian HAtomicSO() const;
    Hamiltonian HRashba(double kx, double ky) const;
};

#endif