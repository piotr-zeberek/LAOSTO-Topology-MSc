#include "ToyModel.h"

void ToyModel::set_nonzero_indices()
{
    std::vector<ElementIndex> all_indices;
    all_indices.reserve(n_bands * n_bands);

    for (int i = 0; i < n_bands; ++i)
    {
        for (int j = 0; j < n_bands; ++j)
        {
            all_indices.emplace_back(i, j);
        }
    }

    std::vector<ElementIndex> Delta_indices = {{0, 1}, {1, 0}};

    Hk_discrete_ky_onsite_nonzero_indices = all_indices;
    Hk_discrete_ky_hopping_p_nonzero_indices = all_indices;
    Hk_discrete_ky_hopping_m_nonzero_indices = all_indices;

    Delta_discrete_ky_onsite_nonzero_indices = Delta_indices;

    Delta_Adjoint_discrete_ky_onsite_nonzero_indices = Delta_indices;

    mHmkT_discrete_ky_onsite_nonzero_indices = all_indices;
    mHmkT_discrete_ky_hopping_p_nonzero_indices = all_indices;
    mHmkT_discrete_ky_hopping_m_nonzero_indices = all_indices;

    Hk_discrete_onsite_nonzero_indices = all_indices;
    Hk_discrete_hopping_xp_nonzero_indices = all_indices;
    Hk_discrete_hopping_xm_nonzero_indices = all_indices;
    Hk_discrete_hopping_yp_nonzero_indices = all_indices;
    Hk_discrete_hopping_ym_nonzero_indices = all_indices;

    Delta_discrete_onsite_nonzero_indices = Delta_indices;

    Delta_Adjoint_discrete_onsite_nonzero_indices = Delta_indices;

    mHmkT_discrete_onsite_nonzero_indices = all_indices;
    mHmkT_discrete_hopping_xp_nonzero_indices = all_indices;
    mHmkT_discrete_hopping_xm_nonzero_indices = all_indices;
    mHmkT_discrete_hopping_yp_nonzero_indices = all_indices;
    mHmkT_discrete_hopping_ym_nonzero_indices = all_indices;
}

Hamiltonian ToyModel::Hk(double kx, double ky) const
{
    return Hkin(kx, ky) + HZeeman() + HRashba(kx, ky) - mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian ToyModel::Delta(double kx, double ky) const
{
    return delta_SC * isy;
}

Hamiltonian ToyModel::Delta_Adjoint(double kx, double ky) const
{
    return (delta_SC * isy).adjoint();
}

Hamiltonian ToyModel::Hk_discrete_ky_onsite(double kx, double y) const
{
    return Hkin(kx) + HZeeman() + HRashba(kx) - mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian ToyModel::Hk_discrete_ky_hopping_p(double kx, double y) const
{
    Hamiltonian H0 = -ty * s0;
    Hamiltonian HRSO = -1i / 2.0 * delta_RSO_y * sx;

    return H0 + HRSO;
}

Hamiltonian ToyModel::Hk_discrete_ky_hopping_m(double kx, double y) const
{
    Hamiltonian H0 = -ty * s0;
    Hamiltonian HRSO = 1i / 2.0 * delta_RSO_y * sx;

    return H0 + HRSO;
}

Hamiltonian ToyModel::Delta_discrete_ky_onsite(double kx, double y) const
{
    return delta_SC * isy;
}

Hamiltonian ToyModel::Delta_discrete_ky_hopping_p(double kx, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian ToyModel::Delta_discrete_ky_hopping_m(double kx, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian ToyModel::Delta_Adjoint_discrete_ky_onsite(double kx, double y) const
{
    return (delta_SC * isy).adjoint();
}

Hamiltonian ToyModel::Delta_Adjoint_discrete_ky_hopping_p(double kx, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian ToyModel::Delta_Adjoint_discrete_ky_hopping_m(double kx, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian ToyModel::mHmkT_discrete_ky_onsite(double kx, double y) const
{
    return -Hkin(kx) - HZeeman().transpose() + HRashba(kx).transpose() + mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian ToyModel::mHmkT_discrete_ky_hopping_p(double kx, double y) const
{
    Hamiltonian H0 = -ty * s0;
    Hamiltonian HRSO = -1i / 2.0 * delta_RSO_y * sx;

    return -H0 + HRSO;
}

Hamiltonian ToyModel::mHmkT_discrete_ky_hopping_m(double kx, double y) const
{
    Hamiltonian H0 = -ty * s0;
    Hamiltonian HRSO = 1i / 2.0 * delta_RSO_y * sx;

    return -H0 + HRSO; // no need to transpose HRSO, symmetric
}

Hamiltonian ToyModel::Hk_discrete_onsite(double x, double y) const
{
    return Hkin() + HZeeman() - mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian ToyModel::Hk_discrete_hopping_xp(double x, double y) const
{
    Hamiltonian H0 = -tx * s0;
    Hamiltonian HRSO = 1i / 2.0 * delta_RSO_x * sy;

    return H0 + HRSO;
}

Hamiltonian ToyModel::Hk_discrete_hopping_xm(double x, double y) const
{
    Hamiltonian H0 = -tx * s0;
    Hamiltonian HRSO = -1i / 2.0 * delta_RSO_x * sy;

    return H0 + HRSO;
}

Hamiltonian ToyModel::Hk_discrete_hopping_yp(double x, double y) const
{
    Hamiltonian H0 = -ty * s0;
    Hamiltonian HRSO = -1i / 2.0 * delta_RSO_y * sx;

    return H0 + HRSO;
}

Hamiltonian ToyModel::Hk_discrete_hopping_ym(double x, double y) const
{
    Hamiltonian H0 = -ty * s0;
    Hamiltonian HRSO = 1i / 2.0 * delta_RSO_y * sx;

    return H0 + HRSO;
}

Hamiltonian ToyModel::Delta_discrete_onsite(double x, double y) const
{
    return delta_SC * isy;
}

Hamiltonian ToyModel::Delta_discrete_hopping_xp(double x, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian ToyModel::Delta_discrete_hopping_xm(double x, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian ToyModel::Delta_discrete_hopping_yp(double x, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian ToyModel::Delta_discrete_hopping_ym(double x, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian ToyModel::Delta_Adjoint_discrete_onsite(double x, double y) const
{
    return (delta_SC * isy).adjoint();
}

Hamiltonian ToyModel::Delta_Adjoint_discrete_hopping_xp(double x, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian ToyModel::Delta_Adjoint_discrete_hopping_xm(double x, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian ToyModel::Delta_Adjoint_discrete_hopping_yp(double x, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian ToyModel::Delta_Adjoint_discrete_hopping_ym(double x, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian ToyModel::mHmkT_discrete_onsite(double x, double y) const
{
    return -Hkin() - HZeeman().transpose() + mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian ToyModel::mHmkT_discrete_hopping_xp(double x, double y) const
{
    return -Hk_discrete_hopping_xp(x, y); // transposing HRSO equivalent to changing sign because of sy
}

Hamiltonian ToyModel::mHmkT_discrete_hopping_xm(double x, double y) const
{
    return -Hk_discrete_hopping_xm(x, y);
}

Hamiltonian ToyModel::mHmkT_discrete_hopping_yp(double x, double y) const
{
    Hamiltonian H0 = -ty * s0;
    Hamiltonian HRSO = -1i / 2.0 * delta_RSO_y * sx;

    return -H0 + HRSO; // no need to transpose HRSO, symmetric
}

Hamiltonian ToyModel::mHmkT_discrete_hopping_ym(double x, double y) const
{
    Hamiltonian H0 = -ty * s0;
    Hamiltonian HRSO = 1i / 2.0 * delta_RSO_y * sx;

    return -H0 + HRSO;
}

Hamiltonian ToyModel::Hkin(double kx, double ky) const
{
    return Ek(kx, ky) * s0;
}

Hamiltonian ToyModel::Hkin(double kx) const
{
    return Ek(kx) * s0;
}

Hamiltonian ToyModel::Hkin() const
{
    return Ek() * s0;
}

Hamiltonian ToyModel::HZeeman() const
{
    Hamiltonian HBx = Bx * sx;
    Hamiltonian HBy = By * sy;
    Hamiltonian HBz = Bz * sz;

    return 0.5 * g_Lande * 0.5 * (HBx + HBy + HBz);
}

Hamiltonian ToyModel::HRashba(double kx, double ky) const
{
    return delta_RSO_y * std::sin(ky) * sx - delta_RSO_x * std::sin(kx) * sy;
}

Hamiltonian ToyModel::HRashba(double kx) const
{
    return -delta_RSO_x * std::sin(kx) * sy;
}
