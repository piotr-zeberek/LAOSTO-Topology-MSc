#ifndef TOYMODEL_H
#define TOYMODEL_H

#include "System2D.h"
#include "utils.h"

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

    ToyModel() : System2D()
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

protected:
    // Triplets only upper triangular part of the entire hamiltonian for onsites
    // continues hamiltonians
    std::vector<Triplet> Hk_triplets(double kx, double ky) const;
    std::vector<Triplet> Delta_triplets(double kx, double ky) const;
    std::vector<Triplet> mHmkT_triplets(double kx, double ky) const;

    // discrete in x direction
    std::vector<Triplet> Hk_discrete_kx_onsite_triplets(double x, double ky) const;
    std::vector<Triplet> Hk_discrete_kx_hopping_p_triplets(double x, double ky) const;

    std::vector<Triplet> Delta_discrete_kx_onsite_triplets(double x, double ky) const;
    // std::vector<Triplet> Delta_discrete_kx_hopping_p_triplets(double x, double ky) const;

    std::vector<Triplet> mHmkT_discrete_kx_onsite_triplets(double x, double ky) const;
    std::vector<Triplet> mHmkT_discrete_kx_hopping_p_triplets(double x, double ky) const;

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
    // std::vector<Triplet> Hk_discrete_hopping_pp_triplets(double x, double y) const;
    // std::vector<Triplet> Hk_discrete_hopping_pm_triplets(double x, double y) const;

    std::vector<Triplet> Delta_discrete_onsite_triplets(double x, double y) const;
    // std::vector<Triplet> Delta_discrete_hopping_xp_triplets(double x, double y) const;
    // std::vector<Triplet> Delta_discrete_hopping_yp_triplets(double x, double y) const;
    // std::vector<Triplet> Delta_discrete_hopping_pp_triplets(double x, double y) const;
    // std::vector<Triplet> Delta_discrete_hopping_pm_triplets(double x, double y) const;

    std::vector<Triplet> mHmkT_discrete_onsite_triplets(double x, double y) const;
    std::vector<Triplet> mHmkT_discrete_hopping_xp_triplets(double x, double y) const;
    std::vector<Triplet> mHmkT_discrete_hopping_yp_triplets(double x, double y) const;
    // std::vector<Triplet> mHmkT_discrete_hopping_pp_triplets(double x, double y) const;
    // std::vector<Triplet> mHmkT_discrete_hopping_pm_triplets(double x, double y) const;

private:
    double Ek(double kx, double ky) const
    {
        return 2.0 * tx * (1.0 - std::cos(kx)) + 2.0 * ty * (1.0 - std::cos(ky));
    }

    double Ek_discrete_kx(double ky) const
    {
        return 2.0 * tx * (1.0 - std::cos(ky)) + 2.0 * ty;
    }

    double Ek_discrete_ky(double kx) const
    {
        return 2.0 * tx * (1.0 - std::cos(kx)) + 2.0 * ty;
    }

    double Ek() const
    {
        return 2.0 * (tx + ty);
    }

    Hamiltonian Hk_mat(double kx, double ky) const;
    Hamiltonian Delta(double kx, double ky) const;

    Hamiltonian Hk_discrete_kx_onsite(double x, double ky) const;
    Hamiltonian Hk_discrete_kx_hopping_p(double x, double ky) const;
    Hamiltonian mHmkT_discrete_kx_onsite(double x, double ky) const;
    Hamiltonian mHmkT_discrete_kx_hopping_p(double x, double ky) const;

    Hamiltonian Hk_discrete_ky_onsite(double kx, double y) const;
    Hamiltonian Hk_discrete_ky_hopping_p(double kx, double y) const;
    Hamiltonian mHmkT_discrete_ky_onsite(double kx, double y) const;
    Hamiltonian mHmkT_discrete_ky_hopping_p(double kx, double y) const;

    Hamiltonian Hk_discrete_onsite(double x, double y) const;
    Hamiltonian Hk_discrete_hopping_xp(double x, double y) const;
    Hamiltonian Hk_discrete_hopping_yp(double x, double y) const;
    Hamiltonian mHmkT_discrete_onsite(double x, double y) const;
    Hamiltonian mHmkT_discrete_hopping_xp(double x, double y) const;
    Hamiltonian mHmkT_discrete_hopping_yp(double x, double y) const;

    Hamiltonian Hkin(double kx, double ky) const;
    Hamiltonian Hkin_discrete_kx(double ky) const;
    Hamiltonian Hkin_discrete_ky(double kx) const;
    Hamiltonian Hkin() const;
    Hamiltonian HZeeman() const;
    Hamiltonian HRashba(double kx, double ky) const;
    Hamiltonian HRashba_discrete_kx(double ky) const;
    Hamiltonian HRashba_discrete_ky(double kx) const;
};

#endif