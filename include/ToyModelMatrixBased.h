#ifndef TOYMODELMATRIXBASED_H
#define TOYMODELMATRIXBASED_H

#include "System2DMatrixBased.h"
#include "utils.h"

struct ToyModelMatrixBased : public System2DMatrixBased
{
    // effective mass
    double m = 1.0;

    // Hopping amplitudes
    double tx = 1.0;
    double ty = 1.0;

    // Rashba coupling
    double delta_RSO_x{};
    double delta_RSO_y{};

    ToyModelMatrixBased()
    {
        set_default_parameters();
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
    Hamiltonian HZeeman(double kx, double ky) const;
    Hamiltonian HRashba(double kx, double ky) const;
};

#endif