#ifndef TOYMODELELEMENTBASED_H
#define TOYMODELELEMENTBASED_H

#include "System2DElementBased.h"

struct ToyModelElementBased : public System2DElementBased
{
    // effective mass
    double m = 1.0;

    // Hopping amplitudes
    double tx = 1.0;
    double ty = 1.0;

    // Rashba coupling
    double delta_RSO_x{};
    double delta_RSO_y{};

    ToyModelElementBased()
    {
        set_default_parameters();

        set_Hk_elements();
        set_Delta_elements();
        set_Delta_Adjoint_elements();
        set_mHmkT_elements();

        set_Hk_discrete_ky_elements();
        set_Delta_discrete_ky_elements();
        set_Delta_Adjoint_discrete_ky_elements();
        set_mHmkT_discrete_ky_elements();

        set_Hk_discrete_elements();
        set_Delta_discrete_elements();
        set_Delta_Adjoint_discrete_elements();
        set_mHmkT_discrete_elements();
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

    void set_Hk_elements() override;
    void set_Delta_elements() override;
    void set_Delta_Adjoint_elements() override;
    void set_mHmkT_elements() override;

    void set_Hk_discrete_ky_elements() override;
    void set_Delta_discrete_ky_elements() override;
    void set_Delta_Adjoint_discrete_ky_elements() override;
    void set_mHmkT_discrete_ky_elements() override;

    void set_Hk_discrete_elements() override;
    void set_Delta_discrete_elements() override;
    void set_Delta_Adjoint_discrete_elements() override;
    void set_mHmkT_discrete_elements() override;

    double Ek(double kx, double ky) const
    {
        return 2.0 * (tx * (1.0 - std::cos(kx)) + ty * (1.0 - std::cos(ky)));
    }

    double Ek(double kx) const
    {
        return 2.0 * (tx * (1.0 - std::cos(kx)) + ty);
    }

    double Ek() const
    {
        return 2.0 * (tx + ty);
    }

    double HBz() const
    {
        return 0.5 * g_Lande * 0.5 * Bz;
    }

    std::complex<double> HBx_iy() const
    {
        return 0.5 * g_Lande * 0.5 * (Bx - 1i * By);
    }

    std::complex<double> HRSO(double kx, double ky) const
    {
        return delta_RSO_y * std::sin(ky) + 1i * delta_RSO_x * std::sin(kx);
    }

    std::complex<double> HRSO(double kx) const
    {
        return 1i * delta_RSO_x * std::sin(kx);
    }
};

#endif