#include "LAOSTO.h"

void LAOSTO::set_nonzero_indices()
{
    for (int i = 0; i < n_bands; ++i)
    {
        for (int j = 0; j < n_bands; ++j)
        {
            if (j != n_bands + 1 - i)
            {
                Hk_discrete_ky_onsite_nonzero_indices.emplace_back(i, j);
            }
        }

        if (i + 2 < n_bands)
        {
            Hk_discrete_ky_hopping_p_nonzero_indices.emplace_back(i, i + 2);
        }

        if (i - 2 >= 0)
        {
            Hk_discrete_ky_hopping_p_nonzero_indices.emplace_back(i, i - 2);
        }

        if (i + 4 < n_bands)
        {
            Hk_discrete_hopping_xp_nonzero_indices.emplace_back(i, i + 4);
        }

        if (i - 4 >= 0)
        {
            Hk_discrete_hopping_xp_nonzero_indices.emplace_back(i, i - 4);
        }

        Hk_discrete_ky_hopping_p_nonzero_indices.emplace_back(i, i);

        Hk_discrete_hopping_xp_nonzero_indices.emplace_back(i, i);

        Hk_discrete_hopping_yp_nonzero_indices.emplace_back(i, i);
    }

    Hk_discrete_hopping_yp_nonzero_indices.emplace_back(0, 2);
    Hk_discrete_hopping_yp_nonzero_indices.emplace_back(1, 3);
    Hk_discrete_hopping_yp_nonzero_indices.emplace_back(2, 0);
    Hk_discrete_hopping_yp_nonzero_indices.emplace_back(3, 1);

    Hk_discrete_onsite_nonzero_indices = Hk_discrete_ky_onsite_nonzero_indices;
    Hk_discrete_ky_hopping_m_nonzero_indices = Hk_discrete_ky_hopping_p_nonzero_indices;
    Hk_discrete_hopping_xm_nonzero_indices = Hk_discrete_hopping_xp_nonzero_indices;
    Hk_discrete_hopping_ym_nonzero_indices = Hk_discrete_hopping_yp_nonzero_indices;

    Delta_discrete_ky_onsite_nonzero_indices.emplace_back(0, 1);
    Delta_discrete_ky_onsite_nonzero_indices.emplace_back(1, 0);
    Delta_discrete_ky_onsite_nonzero_indices.emplace_back(2, 3);
    Delta_discrete_ky_onsite_nonzero_indices.emplace_back(3, 2);
    Delta_discrete_ky_onsite_nonzero_indices.emplace_back(4, 5);
    Delta_discrete_ky_onsite_nonzero_indices.emplace_back(5, 4);

    Delta_discrete_onsite_nonzero_indices = Delta_discrete_ky_onsite_nonzero_indices;

    Delta_Adjoint_discrete_ky_onsite_nonzero_indices = Delta_discrete_ky_onsite_nonzero_indices;
    Delta_Adjoint_discrete_onsite_nonzero_indices = Delta_discrete_onsite_nonzero_indices;

    mHmkT_discrete_ky_onsite_nonzero_indices = Hk_discrete_ky_onsite_nonzero_indices;
    mHmkT_discrete_ky_hopping_p_nonzero_indices = Hk_discrete_ky_hopping_p_nonzero_indices;
    mHmkT_discrete_ky_hopping_m_nonzero_indices = Hk_discrete_ky_hopping_m_nonzero_indices;
    mHmkT_discrete_onsite_nonzero_indices = Hk_discrete_onsite_nonzero_indices;
    mHmkT_discrete_hopping_xp_nonzero_indices = Hk_discrete_hopping_xp_nonzero_indices;
    mHmkT_discrete_hopping_xm_nonzero_indices = Hk_discrete_hopping_xm_nonzero_indices;
    mHmkT_discrete_hopping_yp_nonzero_indices = Hk_discrete_hopping_yp_nonzero_indices;
    mHmkT_discrete_hopping_ym_nonzero_indices = Hk_discrete_hopping_ym_nonzero_indices;
}

Hamiltonian LAOSTO::Hk(double kx, double ky) const
{
    return Hkin(kx, ky) + HZeeman() + HAtomicSO() + HRashba(kx, ky) - mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian LAOSTO::Delta(double kx, double ky) const
{
    return 1i * delta_SC * I3sy;
}

Hamiltonian LAOSTO::Delta_Adjoint(double kx, double ky) const
{
    return (1i * delta_SC * I3sy).adjoint();
}

Hamiltonian LAOSTO::Hk_discrete_ky_onsite(double kx, double y) const
{
    return Hkin(kx) + HZeeman() + HAtomicSO() + HRashba(kx) - mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian LAOSTO::Hk_discrete_ky_hopping_p(double kx, double y) const
{
    Eigen::Matrix3cd ek;
    ek << -tl, 0.0, 0.0,
        0.0, -th, Ek_h(kx),
        0.0, Ek_h(kx), -tl;
    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    Eigen::Matrix3cd rso = Eigen::Matrix3cd::Zero();
    rso(0, 1) = 0.5;
    rso(1, 0) = -0.5;
    Hamiltonian HRSO = delta_RSO * Eigen::kroneckerProduct(rso, s0);

    return H0 + HRSO;
}

Hamiltonian LAOSTO::Hk_discrete_ky_hopping_m(double kx, double y) const
{
    Eigen::Matrix3cd ek;
    ek << -tl, 0.0, 0.0,
        0.0, -th, -Ek_h(kx),
        0.0, -Ek_h(kx), -tl;
    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    Eigen::Matrix3cd rso = Eigen::Matrix3cd::Zero();
    rso(0, 1) = -0.5;
    rso(1, 0) = 0.5;
    Hamiltonian HRSO = delta_RSO * Eigen::kroneckerProduct(rso, s0);

    return H0 + HRSO;
}

Hamiltonian LAOSTO::Delta_discrete_ky_onsite(double kx, double y) const
{
    return 1i * delta_SC * I3sy;
}

Hamiltonian LAOSTO::Delta_discrete_ky_hopping_p(double kx, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian LAOSTO::Delta_discrete_ky_hopping_m(double kx, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian LAOSTO::Delta_Adjoint_discrete_ky_onsite(double kx, double y) const
{
    return (1i * delta_SC * I3sy).adjoint();
}

Hamiltonian LAOSTO::Delta_Adjoint_discrete_ky_hopping_p(double kx, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian LAOSTO::Delta_Adjoint_discrete_ky_hopping_m(double kx, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian LAOSTO::mHmkT_discrete_ky_onsite(double kx, double y) const
{
    return -Hkin(kx) - HZeeman().transpose() - HAtomicSO().transpose() - HRashba(kx) + mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian LAOSTO::mHmkT_discrete_ky_hopping_p(double kx, double y) const
{
    return -Hk_discrete_ky_hopping_p(kx, y);
}

Hamiltonian LAOSTO::mHmkT_discrete_ky_hopping_m(double kx, double y) const
{
    return -Hk_discrete_ky_hopping_m(kx, y);
}

Hamiltonian LAOSTO::Hk_discrete_onsite(double x, double y) const
{
    return Hkin() + HZeeman() + HAtomicSO() - mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian LAOSTO::Hk_discrete_hopping_xp(double x, double y) const
{
    Eigen::Matrix3cd ek;
    ek << -tl, 0.0, 0.0,
        0.0, -tl, 0.0,
        0.0, 0.0, -th;
    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    Eigen::Matrix3cd rso = Eigen::Matrix3cd::Zero();
    rso(0, 2) = 0.5;
    rso(2, 0) = -0.5;
    Hamiltonian HRSO = delta_RSO * Eigen::kroneckerProduct(rso, s0);

    return H0 + HRSO;
}

Hamiltonian LAOSTO::Hk_discrete_hopping_xm(double x, double y) const
{
    Eigen::Matrix3cd ek;
    ek << -tl, 0.0, 0.0,
        0.0, -tl, 0.0,
        0.0, 0.0, -th;
    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    Eigen::Matrix3cd rso = Eigen::Matrix3cd::Zero();
    rso(0, 2) = -0.5;
    rso(2, 0) = 0.5;
    Hamiltonian HRSO = delta_RSO * Eigen::kroneckerProduct(rso, s0);

    return H0 + HRSO;
}

Hamiltonian LAOSTO::Hk_discrete_hopping_yp(double x, double y) const
{
    Eigen::Matrix3cd ek;
    ek << -tl, 0.0, 0.0,
        0.0, -th, 0.0,
        0.0, 0.0, -tl;
    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    Eigen::Matrix3cd rso = Eigen::Matrix3cd::Zero();
    rso(0, 1) = 0.5;
    rso(1, 0) = -0.5;
    Hamiltonian HRSO = delta_RSO * Eigen::kroneckerProduct(rso, s0);

    return H0 + HRSO;
}

Hamiltonian LAOSTO::Hk_discrete_hopping_ym(double x, double y) const
{
    Eigen::Matrix3cd ek;
    ek << -tl, 0.0, 0.0,
        0.0, -th, 0.0,
        0.0, 0.0, -tl;
    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    Eigen::Matrix3cd rso = Eigen::Matrix3cd::Zero();
    rso(0, 1) = -0.5;
    rso(1, 0) = 0.5;
    Hamiltonian HRSO = delta_RSO * Eigen::kroneckerProduct(rso, s0);

    return H0 + HRSO;
}

Hamiltonian LAOSTO::Delta_discrete_onsite(double x, double y) const
{
    return 1i * delta_SC * I3sy;
}

Hamiltonian LAOSTO::Delta_discrete_hopping_xp(double x, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian LAOSTO::Delta_discrete_hopping_xm(double x, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian LAOSTO::Delta_discrete_hopping_yp(double x, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian LAOSTO::Delta_discrete_hopping_ym(double x, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian LAOSTO::Delta_Adjoint_discrete_onsite(double x, double y) const
{
    return (1i * delta_SC * I3sy).adjoint();
}

Hamiltonian LAOSTO::Delta_Adjoint_discrete_hopping_xp(double x, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian LAOSTO::Delta_Adjoint_discrete_hopping_xm(double x, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian LAOSTO::Delta_Adjoint_discrete_hopping_yp(double x, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian LAOSTO::Delta_Adjoint_discrete_hopping_ym(double x, double y) const
{
    return Hamiltonian::Zero(n_bands, n_bands);
}

Hamiltonian LAOSTO::mHmkT_discrete_onsite(double x, double y) const
{
    return -Hkin() - HZeeman().transpose() - HAtomicSO().transpose() + mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian LAOSTO::mHmkT_discrete_hopping_xp(double x, double y) const
{
    return -Hk_discrete_hopping_xp(x, y);
}

Hamiltonian LAOSTO::mHmkT_discrete_hopping_xm(double x, double y) const
{
    return -Hk_discrete_hopping_xm(x, y);
}

Hamiltonian LAOSTO::mHmkT_discrete_hopping_yp(double x, double y) const
{
    return -Hk_discrete_hopping_yp(x, y);
}

Hamiltonian LAOSTO::mHmkT_discrete_hopping_ym(double x, double y) const
{
    return -Hk_discrete_hopping_ym(x, y);
}

Hamiltonian LAOSTO::Hk_discrete_hopping_pp(double x, double y) const
{
    double ek_h = -td / 2.0;

    Eigen::Matrix3cd ek;
    ek << 0.0, 0.0, 0.0,
        0.0, 0.0, ek_h,
        0.0, ek_h, 0.0;

    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    return H0;
}

Hamiltonian LAOSTO::Hk_discrete_hopping_mp(double x, double y) const
{
    return -Hk_discrete_hopping_pp(x, y);
}
Hamiltonian LAOSTO::Hk_discrete_hopping_mm(double x, double y) const
{
    return Hk_discrete_hopping_pp(x, y);
}
Hamiltonian LAOSTO::Hk_discrete_hopping_pm(double x, double y) const
{
    return -Hk_discrete_hopping_pp(x, y);
}

Hamiltonian LAOSTO::mHmkT_discrete_hopping_pp(double x, double y) const
{
    return -Hk_discrete_hopping_pp(x, y);
}

Hamiltonian LAOSTO::mHmkT_discrete_hopping_mp(double x, double y) const
{
    return -Hk_discrete_hopping_mp(x, y);
}

Hamiltonian LAOSTO::mHmkT_discrete_hopping_pm(double x, double y) const
{
    return -Hk_discrete_hopping_pm(x, y);
}

Hamiltonian LAOSTO::mHmkT_discrete_hopping_mm(double x, double y) const
{
    return -Hk_discrete_hopping_mm(x, y);
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

    rso(0, 1) = 1i * std::sin(ky);
    rso(0, 2) = 1i * std::sin(kx);

    rso(1, 0) = -1i * std::sin(ky);
    rso(2, 0) = -1i * std::sin(kx);

    return delta_RSO * Eigen::kroneckerProduct(rso, s0);
}

Hamiltonian LAOSTO::HRashba(double kx) const
{
    Eigen::Matrix3cd rso;

    rso(0, 2) = 1i * std::sin(kx);
    rso(2, 0) = -1i * std::sin(kx);

    return delta_RSO * Eigen::kroneckerProduct(rso, s0);
}

Hamiltonian LAOSTO::HBdG_discrete_hopping_pp(double x, double y) const
{
    Hamiltonian H = Hamiltonian::Zero(2 * n_bands, 2 * n_bands);

    H.block(0, 0, n_bands, n_bands) = Hk_discrete_hopping_pp(x, y);
    H.block(n_bands, n_bands, n_bands, n_bands) = mHmkT_discrete_hopping_pp(x, y);

    return H;
}

Hamiltonian LAOSTO::HBdG_discrete_hopping_mp(double x, double y) const
{
    Hamiltonian H = Hamiltonian::Zero(2 * n_bands, 2 * n_bands);

    H.block(0, 0, n_bands, n_bands) = Hk_discrete_hopping_mp(x, y);
    H.block(n_bands, n_bands, n_bands, n_bands) = mHmkT_discrete_hopping_mp(x, y);

    return H;
}

Hamiltonian LAOSTO::HBdG_discrete_hopping_pm(double x, double y) const
{
    Hamiltonian H = Hamiltonian::Zero(2 * n_bands, 2 * n_bands);

    H.block(0, 0, n_bands, n_bands) = Hk_discrete_hopping_pm(x, y);
    H.block(n_bands, n_bands, n_bands, n_bands) = mHmkT_discrete_hopping_pm(x, y);

    return H;
}

Hamiltonian LAOSTO::HBdG_discrete_hopping_mm(double x, double y) const
{
    Hamiltonian H = Hamiltonian::Zero(2 * n_bands, 2 * n_bands);

    H.block(0, 0, n_bands, n_bands) = Hk_discrete_hopping_mm(x, y);
    H.block(n_bands, n_bands, n_bands, n_bands) = mHmkT_discrete_hopping_mm(x, y);

    return H;
}

Hamiltonian LAOSTO::HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    std::size_t submatrix_size = 2 * n_bands;
    Hamiltonian H = Hamiltonian::Zero(submatrix_size * n_kx * n_ky, submatrix_size * n_kx * n_ky);
    int ix = 0;
    int iy = 0;
    int ixy = 0;

    for (int i = 0; i < n_kx * n_ky; ++i)
    {
        double x = dx * (i / n_ky);
        double y = dy * (i % n_ky);

        Hamiltonian o = HBdG_discrete_onsite(x, y);
        Hamiltonian hxp = HBdG_discrete_hopping_xp(x, y);
        Hamiltonian hxm = HBdG_discrete_hopping_xm(x, y);
        Hamiltonian hyp = HBdG_discrete_hopping_yp(x, y);
        Hamiltonian hym = HBdG_discrete_hopping_ym(x, y);

        Hamiltonian hpp = HBdG_discrete_hopping_pp(x, y);
        Hamiltonian hmp = HBdG_discrete_hopping_mp(x, y);
        Hamiltonian hpm = HBdG_discrete_hopping_pm(x, y);
        Hamiltonian hmm = HBdG_discrete_hopping_mm(x, y);

        H.block(i * submatrix_size, i * submatrix_size, submatrix_size, submatrix_size) = o;

        iy = i - 1;
        if (i % n_ky != 0)
            H.block(i * submatrix_size, iy * submatrix_size, submatrix_size, submatrix_size) = hym;

        iy = i + 1;
        if (iy % n_ky != 0)
            H.block(i * submatrix_size, iy * submatrix_size, submatrix_size, submatrix_size) = hyp;

        ix = i - n_ky;
        if (ix >= 0)
            H.block(i * submatrix_size, ix * submatrix_size, submatrix_size, submatrix_size) = hxm;

        ix = i + n_ky;
        if (ix < n_kx * n_ky)
            H.block(i * submatrix_size, ix * submatrix_size, submatrix_size, submatrix_size) = hxp;

        ixy = i - n_ky - 1;
        if (ixy >= 0 && i % n_ky != 0)
            H.block(i * submatrix_size, ixy * submatrix_size, submatrix_size, submatrix_size) = hmm;

        ixy = i - n_ky + 1;
        if (ixy >= 0 && ixy % n_ky != 0)
            H.block(i * submatrix_size, ixy * submatrix_size, submatrix_size, submatrix_size) = hmp;

        ixy = i + n_ky - 1;
        if (ixy < n_kx * n_ky && i % n_ky != 0)
            H.block(i * submatrix_size, ixy * submatrix_size, submatrix_size, submatrix_size) = hpm;

        ixy = i + n_ky + 1;
        if (ixy < n_kx * n_ky && ixy % n_ky != 0)
            H.block(i * submatrix_size, ixy * submatrix_size, submatrix_size, submatrix_size) = hpp;
    }

    return H;
}

std::vector<Triplet> LAOSTO::triplets_HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    std::size_t submatrix_size = 2 * n_bands;
    std::size_t n_elements = HBdG_discrete_onsite_nonzero_indices.size() +
                             HBdG_discrete_hopping_xp_nonzero_indices.size() +
                             HBdG_discrete_hopping_xm_nonzero_indices.size() +
                             HBdG_discrete_hopping_yp_nonzero_indices.size() +
                             HBdG_discrete_hopping_ym_nonzero_indices.size() +
                             HBdG_discrete_hopping_pp_nonzero_indices.size() +
                             HBdG_discrete_hopping_mp_nonzero_indices.size() +
                             HBdG_discrete_hopping_pm_nonzero_indices.size() +
                             HBdG_discrete_hopping_mm_nonzero_indices.size();
    std::vector<Triplet> triplets;
    triplets.reserve(n_elements * n_kx * n_ky);
    int ix = 0;
    int iy = 0;
    int ixy = 0;

    for (int i = 0; i < n_kx * n_ky; ++i)
    {
        double x = dx * (i / n_ky);
        double y = dy * (i % n_ky);

        for (const auto &index : HBdG_discrete_onsite_nonzero_indices)
        {
            auto [row, col] = index;
            triplets.emplace_back(i * submatrix_size + row, i * submatrix_size + col, HBdG_discrete_onsite(x, y)(row, col));
        }

        iy = i - 1;
        if (i % n_ky != 0)
            for (const auto &index : HBdG_discrete_hopping_ym_nonzero_indices)
            {
                auto [row, col] = index;
                triplets.emplace_back(i * submatrix_size + row, iy * submatrix_size + col, HBdG_discrete_hopping_ym(x, y)(row, col));
            }

        iy = i + 1;
        if (iy % n_ky != 0)
            for (const auto &index : HBdG_discrete_hopping_yp_nonzero_indices)
            {
                auto [row, col] = index;
                triplets.emplace_back(i * submatrix_size + row, iy * submatrix_size + col, HBdG_discrete_hopping_yp(x, y)(row, col));
            }

        ix = i - n_ky;
        if (ix >= 0)
            for (const auto &index : HBdG_discrete_hopping_xm_nonzero_indices)
            {
                auto [row, col] = index;
                triplets.emplace_back(i * submatrix_size + row, ix * submatrix_size + col, HBdG_discrete_hopping_xm(x, y)(row, col));
            }

        ix = i + n_ky;
        if (ix < n_kx * n_ky)
            for (const auto &index : HBdG_discrete_hopping_xp_nonzero_indices)
            {
                auto [row, col] = index;
                triplets.emplace_back(i * submatrix_size + row, ix * submatrix_size + col, HBdG_discrete_hopping_xp(x, y)(row, col));
            }

        ixy = i - n_ky - 1;
        if (ixy >= 0 && i % n_ky != 0)
            for (const auto &index : HBdG_discrete_hopping_mm_nonzero_indices)
            {
                auto [row, col] = index;
                triplets.emplace_back(i * submatrix_size + row, ixy * submatrix_size + col, HBdG_discrete_hopping_mm(x, y)(row, col));
            }

        ixy = i - n_ky + 1;
        if (ixy >= 0 && ixy % n_ky != 0)
            for (const auto &index : HBdG_discrete_hopping_mp_nonzero_indices)
            {
                auto [row, col] = index;
                triplets.emplace_back(i * submatrix_size + row, ixy * submatrix_size + col, HBdG_discrete_hopping_mp(x, y)(row, col));
            }

        ixy = i + n_ky - 1;
        if (ixy < n_kx * n_ky && i % n_ky != 0)
            for (const auto &index : HBdG_discrete_hopping_pm_nonzero_indices)
            {
                auto [row, col] = index;
                triplets.emplace_back(i * submatrix_size + row, ixy * submatrix_size + col, HBdG_discrete_hopping_pm(x, y)(row, col));
            }

        ixy = i + n_ky + 1;
        if (ixy < n_kx * n_ky && ixy % n_ky != 0)
            for (const auto &index : HBdG_discrete_hopping_pp_nonzero_indices)
            {
                auto [row, col] = index;
                triplets.emplace_back(i * submatrix_size + row, ixy * submatrix_size + col, HBdG_discrete_hopping_pp(x, y)(row, col));
            }
    }

    return triplets;
}