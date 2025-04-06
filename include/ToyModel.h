#ifndef TOYMODEL_H
#define TOYMODEL_H

#include "System2D.h"
#include "utils.h"

#include <unsupported/Eigen/KroneckerProduct>

struct ToyModel : public System2D
{
    // effective mass
    double m = 1.0;

    // Hopping amplitudes
    double tx = 1.0;
    double ty = 1.0;

    // Rashba coupling
    double delta_RSO_x{};
    double delta_RSO_y{};

    // for Delta
    Hamiltonian isy = 1i * sy;

    ToyModel()
    {
        set_default_parameters();
        set_nonzero_indices();
        update_HBdG_nonzero_indices();
    }

    void set_default_parameters() override
    {
        n_bands = 2;

        double dx_nm = 8.0;
        double dy_nm = 8.0;

        m = 0.014;
        dx = nm2au(dx_nm);
        dy = nm2au(dy_nm);
        update_kinetic_hopping_amplitudes();

        double alpha = 50.0; // Rashba coupling in meV nm
        delta_RSO_x = meV2au(alpha / dx_nm);
        delta_RSO_y = meV2au(alpha / dy_nm);

        mu = 0.0;

        g_Lande = -50;
        Bx = 0.0;
        By = 0.0;
        Bz = 0.0;

        delta_SC = meV2au(0.2);
    }

    void update_kinetic_hopping_amplitudes()
    {
        tx = 1.0 / (2.0 * m * dx * dx);
        ty = 1.0 / (2.0 * m * dy * dy);
    }

    void set_nonzero_indices() override;

    Hamiltonian Hk(double kx, double ky) const override;
    Hamiltonian Delta(double kx, double ky) const override;
    Hamiltonian Delta_Adjoint(double kx, double ky) const override;

    // hamiltonian elements

    // discretized in ky
    Hamiltonian Hk_discrete_ky_onsite(double kx, double y) const override;
    Hamiltonian Hk_discrete_ky_hopping_p(double kx, double y) const override;
    Hamiltonian Hk_discrete_ky_hopping_m(double kx, double y) const override;

    Hamiltonian Delta_discrete_ky_onsite(double kx, double y) const override;
    Hamiltonian Delta_discrete_ky_hopping_p(double kx, double y) const override;
    Hamiltonian Delta_discrete_ky_hopping_m(double kx, double y) const override;

    Hamiltonian Delta_Adjoint_discrete_ky_onsite(double kx, double y) const override;
    Hamiltonian Delta_Adjoint_discrete_ky_hopping_p(double kx, double y) const override;
    Hamiltonian Delta_Adjoint_discrete_ky_hopping_m(double kx, double y) const override;

    Hamiltonian mHmkT_discrete_ky_onsite(double kx, double y) const override;
    Hamiltonian mHmkT_discrete_ky_hopping_p(double kx, double y) const override;
    Hamiltonian mHmkT_discrete_ky_hopping_m(double kx, double y) const override;

    // discretized in kx,ky
    Hamiltonian Hk_discrete_onsite(double x, double y) const override;
    Hamiltonian Hk_discrete_hopping_xp(double x, double y) const override;
    Hamiltonian Hk_discrete_hopping_xm(double x, double y) const override;
    Hamiltonian Hk_discrete_hopping_yp(double x, double y) const override;
    Hamiltonian Hk_discrete_hopping_ym(double x, double y) const override;

    Hamiltonian Delta_discrete_onsite(double x, double y) const override;
    Hamiltonian Delta_discrete_hopping_xp(double x, double y) const override;
    Hamiltonian Delta_discrete_hopping_xm(double x, double y) const override;
    Hamiltonian Delta_discrete_hopping_yp(double x, double y) const override;
    Hamiltonian Delta_discrete_hopping_ym(double x, double y) const override;

    Hamiltonian Delta_Adjoint_discrete_onsite(double x, double y) const override;
    Hamiltonian Delta_Adjoint_discrete_hopping_xp(double x, double y) const override;
    Hamiltonian Delta_Adjoint_discrete_hopping_xm(double x, double y) const override;
    Hamiltonian Delta_Adjoint_discrete_hopping_yp(double x, double y) const override;
    Hamiltonian Delta_Adjoint_discrete_hopping_ym(double x, double y) const override;

    Hamiltonian mHmkT_discrete_onsite(double x, double y) const override;
    Hamiltonian mHmkT_discrete_hopping_xp(double x, double y) const override;
    Hamiltonian mHmkT_discrete_hopping_xm(double x, double y) const override;
    Hamiltonian mHmkT_discrete_hopping_yp(double x, double y) const override;
    Hamiltonian mHmkT_discrete_hopping_ym(double x, double y) const override;

    Hamiltonian Hkin(double kx, double ky) const;
    Hamiltonian Hkin(double kx) const;
    Hamiltonian Hkin() const;
    Hamiltonian HZeeman() const;
    Hamiltonian HRashba(double kx, double ky) const;
    Hamiltonian HRashba(double kx) const;

    double Ek(double kx, double ky) const
    {
        return 2.0 * tx * (1.0 - std::cos(kx)) + 2.0 * ty * (1.0 - std::cos(ky));
    }

    double Ek(double kx) const
    {
        return 2.0 * tx * (1.0 - std::cos(kx)) + 2.0 * ty;
    }

    double Ek() const
    {
        return 2.0 * (tx + ty);
    }
};

#endif