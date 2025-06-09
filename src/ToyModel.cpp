#include "ToyModel.h"

std::vector<Triplet> ToyModel::Hk_triplets(double kx, double ky) const
{
    return generate_triplets(Hk_mat(kx, ky), true);
}

std::vector<Triplet> ToyModel::Delta_triplets(double kx, double ky) const
{
    return {
        {0, 1, delta_SC},
        {1, 0, -delta_SC}};
}

std::vector<Triplet> ToyModel::mHmkT_triplets(double kx, double ky) const
{
    return generate_triplets(-Hk_mat(-kx, -ky).transpose(), true);
}

std::vector<Triplet> ToyModel::Hk_discrete_ky_onsite_triplets(double kx, double y) const
{
    return generate_triplets(Hk_discrete_ky_onsite(kx, y), true);
}

std::vector<Triplet> ToyModel::Hk_discrete_ky_hopping_p_triplets(double kx, double y) const
{
    return generate_triplets(Hk_discrete_ky_hopping_p(kx, y));
}

std::vector<Triplet> ToyModel::Delta_discrete_ky_onsite_triplets(double kx, double y) const
{
    return Delta_triplets(kx, y);
}

std::vector<Triplet> ToyModel::mHmkT_discrete_ky_onsite_triplets(double kx, double y) const
{
    return generate_triplets(mHmkT_discrete_ky_onsite(kx, y), true);
}

std::vector<Triplet> ToyModel::mHmkT_discrete_ky_hopping_p_triplets(double kx, double y) const
{
    return generate_triplets(mHmkT_discrete_ky_hopping_p(kx, y));
}

std::vector<Triplet> ToyModel::Hk_discrete_onsite_triplets(double x, double y) const
{
    return generate_triplets(Hk_discrete_onsite(x, y), true);
}

std::vector<Triplet> ToyModel::Hk_discrete_hopping_xp_triplets(double x, double y) const
{
    return generate_triplets(Hk_discrete_hopping_xp(x, y));
}

std::vector<Triplet> ToyModel::Hk_discrete_hopping_yp_triplets(double x, double y) const
{
    return generate_triplets(Hk_discrete_hopping_yp(x, y));
}

std::vector<Triplet> ToyModel::Delta_discrete_onsite_triplets(double x, double y) const
{
    return Delta_triplets(x, y);
}

std::vector<Triplet> ToyModel::mHmkT_discrete_onsite_triplets(double x, double y) const
{
    return generate_triplets(mHmkT_discrete_onsite(x, y), true);
}

std::vector<Triplet> ToyModel::mHmkT_discrete_hopping_xp_triplets(double x, double y) const
{
    return generate_triplets(mHmkT_discrete_hopping_xp(x, y));
}

std::vector<Triplet> ToyModel::mHmkT_discrete_hopping_yp_triplets(double x, double y) const
{
    return generate_triplets(mHmkT_discrete_hopping_yp(x, y));
}

Hamiltonian ToyModel::Hk_mat(double kx, double ky) const
{
    return Hkin(kx, ky) + HZeeman() + HRashba(kx, ky) - mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian ToyModel::Delta(double kx, double ky) const
{
    Hamiltonian isy = 1i * sy;
    return delta_SC * isy;
}

Hamiltonian ToyModel::Hk_discrete_ky_onsite(double kx, double y) const
{
    return Hkin(kx) + HZeeman() + HRashba(kx) - mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian ToyModel::Hk_discrete_ky_hopping_p(double kx, double y) const
{
    Hamiltonian H0 = -ty * s0;
    Hamiltonian HRSO = 0.5 * 1i * delta_RSO_y * sx;

    return H0 + HRSO;
}

Hamiltonian ToyModel::mHmkT_discrete_ky_onsite(double kx, double y) const
{
    return -Hkin(kx) - HZeeman().transpose() + HRashba(kx).transpose() + mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian ToyModel::mHmkT_discrete_ky_hopping_p(double kx, double y) const
{
    Hamiltonian H0 = -ty * s0;
    Hamiltonian HRSO = 0.5 * 1i * delta_RSO_y * sx;

    return -H0 + HRSO;
}

Hamiltonian ToyModel::Hk_discrete_onsite(double x, double y) const
{
    return Hkin() + HZeeman() - mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian ToyModel::Hk_discrete_hopping_xp(double x, double y) const
{
    Hamiltonian H0 = -tx * s0;
    Hamiltonian HRSO = -0.5 * 1i * delta_RSO_x * sy;

    return H0 + HRSO;
}

Hamiltonian ToyModel::Hk_discrete_hopping_yp(double x, double y) const
{
    return Hk_discrete_ky_hopping_p(x, y); 
}

Hamiltonian ToyModel::mHmkT_discrete_onsite(double x, double y) const
{
    return -Hkin() - HZeeman().transpose() + mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian ToyModel::mHmkT_discrete_hopping_xp(double x, double y) const
{
    return -Hk_discrete_hopping_xp(x, y); // transposing HRSO equivalent to changing sign because of sy
}

Hamiltonian ToyModel::mHmkT_discrete_hopping_yp(double x, double y) const
{
    Hamiltonian H0 = -ty * s0;
    Hamiltonian HRSO = 0.5 * 1i * delta_RSO_y * sx;

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
