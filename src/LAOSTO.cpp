#include "LAOSTO.h"

std::vector<Triplet> LAOSTO::Hk_triplets(double kx, double ky) const
{
    return generate_triplets(Hk_mat(kx, ky), true);
}

std::vector<Triplet> LAOSTO::Delta_triplets(double kx, double ky) const
{
    return {
        {0, 1, delta_SC},
        {1, 0, -delta_SC},
        {2, 3, delta_SC},
        {3, 2, -delta_SC},
        {4, 5, delta_SC},
        {5, 4, -delta_SC}};
}

std::vector<Triplet> LAOSTO::mHmkT_triplets(double kx, double ky) const
{
    return generate_triplets(-Hk_mat(-kx, -ky).transpose(), true);
}

std::vector<Triplet> LAOSTO::Hk_discrete_ky_onsite_triplets(double kx, double y) const
{
    return generate_triplets(Hk_discrete_ky_onsite(kx, y), true);
}

std::vector<Triplet> LAOSTO::Hk_discrete_ky_hopping_p_triplets(double kx, double y) const
{
    return {
        // kin
        {0, 0, -tl},
        {1, 1, -tl},
        {2, 2, -th},
        {3, 3, -th},
        {4, 4, -tl},
        {5, 5, -tl},
        {2, 4, Ek_h(kx)},
        {4, 2, Ek_h(kx)},
        {3, 5, Ek_h(kx)},
        {5, 3, Ek_h(kx)},
        // rso
        {0, 2, -0.5 * delta_RSO},
        {2, 0, 0.5 * delta_RSO},
        {1, 3, -0.5 * delta_RSO},
        {3, 1, 0.5 * delta_RSO}};
}

std::vector<Triplet> LAOSTO::Delta_discrete_ky_onsite_triplets(double kx, double y) const
{
    return Delta_triplets(kx, y);
}

std::vector<Triplet> LAOSTO::mHmkT_discrete_ky_onsite_triplets(double kx, double y) const
{
    return generate_triplets(mHmkT_discrete_ky_onsite(kx, y), true);
}

std::vector<Triplet> LAOSTO::mHmkT_discrete_ky_hopping_p_triplets(double kx, double y) const
{
    return negate_triplets(Hk_discrete_ky_hopping_p_triplets(kx, y));
}

std::vector<Triplet> LAOSTO::Hk_discrete_onsite_triplets(double x, double y) const
{
    return generate_triplets(Hk_discrete_onsite(x, y), true);
}

std::vector<Triplet> LAOSTO::Hk_discrete_hopping_xp_triplets(double x, double y) const
{
    return {
        // kin
        {0, 0, -tl},
        {1, 1, -tl},
        {2, 2, -tl},
        {3, 3, -tl},
        {4, 4, -th},
        {5, 5, -th},
        // rso
        {0, 4, -0.5 * delta_RSO},
        {4, 0, 0.5 * delta_RSO},
        {1, 5, -0.5 * delta_RSO},
        {5, 1, 0.5 * delta_RSO}};
}

std::vector<Triplet> LAOSTO::Hk_discrete_hopping_yp_triplets(double x, double y) const
{
    return {
        // kin
        {0, 0, -tl},
        {1, 1, -tl},
        {2, 2, -th},
        {3, 3, -th},
        {4, 4, -tl},
        {5, 5, -tl},
        // rso
        {0, 2, -0.5 * delta_RSO},
        {2, 0, 0.5 * delta_RSO},
        {1, 3, -0.5 * delta_RSO},
        {3, 1, 0.5 * delta_RSO}};
}

std::vector<Triplet> LAOSTO::Hk_discrete_hopping_pp_triplets(double x, double y) const
{
    return {
        // kin
        {2, 4, -0.5 * td},
        {4, 2, -0.5 * td},
        {3, 5, -0.5 * td},
        {5, 3, -0.5 * td}};
}

std::vector<Triplet> LAOSTO::Hk_discrete_hopping_pm_triplets(double x, double y) const
{
    return {
        // kin
        {2, 4, 0.5 * td},
        {4, 2, 0.5 * td},
        {3, 5, 0.5 * td},
        {5, 3, 0.5 * td}};
}

std::vector<Triplet> LAOSTO::Delta_discrete_onsite_triplets(double x, double y) const {
    return Delta_triplets(x, y);
}

std::vector<Triplet> LAOSTO::mHmkT_discrete_onsite_triplets(double x, double y) const {
    return generate_triplets(mHmkT_discrete_onsite(x, y), true);
}

std::vector<Triplet> LAOSTO::mHmkT_discrete_hopping_xp_triplets(double x, double y) const {
    return negate_triplets(Hk_discrete_hopping_xp_triplets(x, y));
}

std::vector<Triplet> LAOSTO::mHmkT_discrete_hopping_yp_triplets(double x, double y) const {
    return negate_triplets(Hk_discrete_hopping_yp_triplets(x, y));
}

std::vector<Triplet> LAOSTO::mHmkT_discrete_hopping_pp_triplets(double x, double y) const {
    return negate_triplets(Hk_discrete_hopping_pp_triplets(x, y));
}

std::vector<Triplet> LAOSTO::mHmkT_discrete_hopping_pm_triplets(double x, double y) const {
    return negate_triplets(Hk_discrete_hopping_pm_triplets(x, y));
}

Hamiltonian LAOSTO::Hk_mat(double kx, double ky) const
{
    return Hkin(kx, ky) + HZeeman() + HAtomicSO() + HRashba(kx, ky) - mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian LAOSTO::Delta(double kx, double ky) const
{
    return 1i * delta_SC * I3sy;
}

Hamiltonian LAOSTO::Hk_discrete_ky_onsite(double kx, double y) const
{
    return Hkin(kx) + HZeeman() + HAtomicSO() + HRashba(kx) - mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian LAOSTO::mHmkT_discrete_ky_onsite(double kx, double y) const
{
    return -Hkin(kx) - HZeeman().transpose() - HAtomicSO().transpose() - HRashba(kx) + mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian LAOSTO::Hk_discrete_onsite(double x, double y) const
{
    return Hkin() + HZeeman() + HAtomicSO() - mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian LAOSTO::mHmkT_discrete_onsite(double x, double y) const
{
    return -Hkin() - HZeeman().transpose() - HAtomicSO().transpose() + mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian LAOSTO::Hkin(double kx, double ky) const
{
    Eigen::Matrix3cd ek;
    ek << Ekxy(kx, ky), 0.0, 0.0,
        0.0, Ekxz(kx, ky), Ek_h(kx, ky),
        0.0, Ek_h(kx, ky), Ekyz(kx, ky);

    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    return Eigen::kroneckerProduct(ek, s0);
}

Hamiltonian LAOSTO::Hkin(double kx) const
{
    Eigen::Matrix3cd ek;
    ek << Ekxy(kx), 0.0, 0.0,
        0.0, Ekxz(kx), 0.0,
        0.0, 0.0, Ekyz(kx);

    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    return H0;
}

Hamiltonian LAOSTO::Hkin() const
{
    Eigen::Matrix3cd ek;
    ek << Ekxy(), 0.0, 0.0,
        0.0, Ekxz(), 0.0,
        0.0, 0.0, Ekyz();

    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    return H0;
}

Hamiltonian LAOSTO::HZeeman() const
{
    Hamiltonian HBx = Bx * (Lxs0 + 0.5 * g_Lande * I3sx);
    Hamiltonian HBy = By * (Lys0 + 0.5 * g_Lande * I3sy);
    Hamiltonian HBz = Bz * (Lzs0 + 0.5 * g_Lande * I3sz);

    return 0.5 * (HBx + HBy + HBz);
}

Hamiltonian LAOSTO::HAtomicSO() const
{
    return delta_SO / 3.0 * HSO_mat;
}

Hamiltonian LAOSTO::HRashba(double kx, double ky) const
{
    Eigen::Matrix3cd rso;
    rso.setZero();

    rso(0, 1) = 1i * std::sin(ky);
    rso(0, 2) = 1i * std::sin(kx);

    rso(1, 0) = -1i * std::sin(ky);
    rso(2, 0) = -1i * std::sin(kx);

    return delta_RSO * Eigen::kroneckerProduct(rso, s0);
}

Hamiltonian LAOSTO::HRashba(double kx) const
{
    Eigen::Matrix3cd rso;
    rso.setZero();

    rso(0, 2) = 1i * std::sin(kx);
    rso(2, 0) = -1i * std::sin(kx);

    return delta_RSO * Eigen::kroneckerProduct(rso, s0);
}
