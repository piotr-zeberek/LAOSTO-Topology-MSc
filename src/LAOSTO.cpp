#include "LAOSTO.h"


Hamiltonian LAOSTO::Hk(double kx, double ky) const
{
    return Hkin(kx, ky) + HZeeman() + HAtomicSO() + HRashba(kx, ky) - mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
}

Hamiltonian LAOSTO::HBdG(double kx, double ky) const
{
    Eigen::MatrixXcd res = -delta_SC * syI3sy;

    res.topLeftCorner(n_bands, n_bands) = Hk(kx, ky);
    res.bottomRightCorner(n_bands, n_bands) = -Hk(-kx, -ky).transpose();

    return res;
}

Hamiltonian LAOSTO::HBdG_discrete_ky_onsite(double kx, double y) const
{
    double ek_xy = 2 * tl * (2 - std::cos(kx)) - delta_E;
    double ek_xz = 2 * tl * (1 - std::cos(kx)) + 2 * th;
    double ek_yz = 2 * tl + 2 * th * (1 - std::cos(kx));

    Eigen::Matrix3cd ek;
    ek << ek_xy, 0.0, 0.0,
        0.0, ek_xz, 0.0,
        0.0, 0.0, ek_yz;

    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    Eigen::Matrix3cd rso = Eigen::Matrix3cd::Zero();
    rso(0, 2) = 1i * std::sin(kx);
    rso(2, 0) = -1i * std::sin(kx);
    Hamiltonian HRSO = delta_RSO * Eigen::kroneckerProduct(rso, s0);

    Hamiltonian HSO = HAtomicSO();
    Hamiltonian HB = HZeeman();

    Hamiltonian res = -delta_SC * syI3sy;

    res.topLeftCorner(n_bands, n_bands) = H0 + HSO + HRSO + HB - mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
    res.bottomRightCorner(n_bands, n_bands) = -H0 - HSO.transpose() - HRSO - HB.transpose() + mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);

    return res;
}

Hamiltonian LAOSTO::HBdG_discrete_ky_hopping_p(double kx, double y) const
{
    double ek_xy = -tl;
    double ek_xz = -th;
    double ek_yz = -tl;

    std::complex<double> ek_h = -1i * td * std::sin(kx);

    Eigen::Matrix3cd ek;
    ek << ek_xy, 0.0, 0.0,
        0.0, ek_xz, ek_h,
        0.0, ek_h, ek_yz;

    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    // kinetic -H(-k)^T
    Eigen::Matrix3cd ek_br; // bottom right
    ek_br << ek_xy, 0.0, 0.0,
        0.0, ek_xz, -ek_h,
        0.0, -ek_h, ek_yz; // minus ek_h, because of -k

    Hamiltonian H0_br = Eigen::kroneckerProduct(ek_br, s0); // no need to transpose,

    Eigen::Matrix3cd rso = Eigen::Matrix3cd::Zero();
    rso(0, 1) = 0.5;
    rso(1, 0) = -0.5;
    Hamiltonian HRSO = delta_RSO * Eigen::kroneckerProduct(rso, s0);

    Hamiltonian res = Hamiltonian::Zero(n_bands_sc, n_bands_sc);

    res.topLeftCorner(n_bands, n_bands) = H0 + HRSO;
    res.bottomRightCorner(n_bands, n_bands) = -H0_br + HRSO; // HRSO transpose and minus cancels out

    return res;
}

Hamiltonian LAOSTO::HBdG_discrete_ky_hopping_m(double kx, double y) const
{
    double ek_xy = -tl;
    double ek_xz = -th;
    double ek_yz = -tl;

    std::complex<double> ek_h = 1i * td * std::sin(kx);

    Eigen::Matrix3cd ek;
    ek << ek_xy, 0.0, 0.0,
        0.0, ek_xz, ek_h,
        0.0, ek_h, ek_yz;

    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    // kinetic -H(-k)^T
    Eigen::Matrix3cd ek_br; // bottom right
    ek_br << ek_xy, 0.0, 0.0,
        0.0, ek_xz, -ek_h,
        0.0, -ek_h, ek_yz; // minus ek_h, because of -k

    Hamiltonian H0_br = Eigen::kroneckerProduct(ek_br, s0); // no need to transpose,

    Eigen::Matrix3cd rso = Eigen::Matrix3cd::Zero();
    rso(0, 1) = -0.5;
    rso(1, 0) = 0.5;
    Hamiltonian HRSO = delta_RSO * Eigen::kroneckerProduct(rso, s0);

    Hamiltonian res = Hamiltonian::Zero(n_bands_sc, n_bands_sc);

    res.topLeftCorner(n_bands, n_bands) = H0 + HRSO;
    res.bottomRightCorner(n_bands, n_bands) = -H0_br + HRSO; // HRSO transpose and minus cancels out

    return res;
}

Hamiltonian LAOSTO::HBdG_discrete_onsite(double x, double y) const
{
    double ek_xy = 4 * tl - delta_E;
    double ek_xz = 2 * (tl + th);
    double ek_yz = 2 * (tl + th);

    Eigen::Matrix3cd ek;
    ek << ek_xy, 0.0, 0.0,
        0.0, ek_xz, 0.0,
        0.0, 0.0, ek_yz;

    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    Hamiltonian HSO = HAtomicSO();
    Hamiltonian HB = HZeeman();

    Hamiltonian res = -delta_SC * syI3sy;

    res.topLeftCorner(n_bands, n_bands) = H0 + HSO + HB - mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);
    res.bottomRightCorner(n_bands, n_bands) = -H0 - HSO.transpose() - HB.transpose() + mu * Eigen::MatrixXcd::Identity(n_bands, n_bands);

    return res;
}

Hamiltonian LAOSTO::HBdG_discrete_hopping_xp(double x, double y) const
{
    double ek_xy = -tl;
    double ek_xz = -tl; // tl for hopping_x, th for hopping_y
    double ek_yz = -th; // th for hopping_x, tl for hopping_y

    Eigen::Matrix3cd ek;
    ek << ek_xy, 0.0, 0.0,
        0.0, ek_xz, 0.0,
        0.0, 0.0, ek_yz;

    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    Eigen::Matrix3cd rso = Eigen::Matrix3cd::Zero();
    rso(0, 2) = 0.5;
    rso(2, 0) = -0.5;
    Hamiltonian HRSO = delta_RSO * Eigen::kroneckerProduct(rso, s0);

    Hamiltonian res = Hamiltonian::Zero(n_bands_sc, n_bands_sc);

    res.topLeftCorner(n_bands, n_bands) = H0 + HRSO;
    res.bottomRightCorner(n_bands, n_bands) = -H0 + HRSO; // HRSO transpose and minus cancels out

    return res;
}

Hamiltonian LAOSTO::HBdG_discrete_hopping_xm(double x, double y) const
{
    double ek_xy = -tl;
    double ek_xz = -tl; // tl for hopping_x, th for hopping_y
    double ek_yz = -th; // th for hopping_x, tl for hopping_y

    Eigen::Matrix3cd ek;
    ek << ek_xy, 0.0, 0.0,
        0.0, ek_xz, 0.0,
        0.0, 0.0, ek_yz;

    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    Eigen::Matrix3cd rso = Eigen::Matrix3cd::Zero();
    rso(0, 2) = -0.5;
    rso(2, 0) = 0.5;
    Hamiltonian HRSO = delta_RSO * Eigen::kroneckerProduct(rso, s0);

    Hamiltonian res = Hamiltonian::Zero(n_bands_sc, n_bands_sc);

    res.topLeftCorner(n_bands, n_bands) = H0 + HRSO;
    res.bottomRightCorner(n_bands, n_bands) = -H0 + HRSO; // HRSO transpose and minus cancels out

    return res;
}

Hamiltonian LAOSTO::HBdG_discrete_hopping_yp(double x, double y) const
{
    double ek_xy = -tl;
    double ek_xz = -th; // tl for hopping_x, th for hopping_y
    double ek_yz = -tl; // th for hopping_x, tl for hopping_y

    Eigen::Matrix3cd ek;
    ek << ek_xy, 0.0, 0.0,
        0.0, ek_xz, 0.0,
        0.0, 0.0, ek_yz;

    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    Eigen::Matrix3cd rso = Eigen::Matrix3cd::Zero();
    rso(0, 1) = 0.5;
    rso(1, 0) = -0.5;
    Hamiltonian HRSO = delta_RSO * Eigen::kroneckerProduct(rso, s0);

    Hamiltonian res = Hamiltonian::Zero(n_bands_sc, n_bands_sc);

    res.topLeftCorner(n_bands, n_bands) = H0 + HRSO;
    res.bottomRightCorner(n_bands, n_bands) = -H0 + HRSO; // HRSO transpose and minus cancels out

    return res;
}

Hamiltonian LAOSTO::HBdG_discrete_hopping_ym(double x, double y) const
{
    double ek_xy = -tl;
    double ek_xz = -th; // tl for hopping_x, th for hopping_y
    double ek_yz = -tl; // th for hopping_x, tl for hopping_y

    Eigen::Matrix3cd ek;
    ek << ek_xy, 0.0, 0.0,
        0.0, ek_xz, 0.0,
        0.0, 0.0, ek_yz;

    Hamiltonian H0 = Eigen::kroneckerProduct(ek, s0);

    Eigen::Matrix3cd rso = Eigen::Matrix3cd::Zero();
    rso(0, 1) = -0.5;
    rso(1, 0) = 0.5;
    Hamiltonian HRSO = delta_RSO * Eigen::kroneckerProduct(rso, s0);

    Hamiltonian res = Hamiltonian::Zero(n_bands_sc, n_bands_sc);

    res.topLeftCorner(n_bands, n_bands) = H0 + HRSO;
    res.bottomRightCorner(n_bands, n_bands) = -H0 + HRSO; // HRSO transpose and minus cancels out

    return res;
}

Hamiltonian LAOSTO::Hkin(double kx, double ky) const
{
    double ek_xy = 2 * tl * (2 - std::cos(kx) - std::cos(ky)) - delta_E;
    double ek_xz = 2 * tl * (1 - std::cos(kx)) + 2 * th * (1 - std::cos(ky));
    double ek_yz = 2 * tl * (1 - std::cos(ky)) + 2 * th * (1 - std::cos(kx));
    double ek_h = 2 * td * std::sin(kx) * std::sin(ky);

    Eigen::Matrix3cd ek;
    ek << ek_xy, 0.0, 0.0,
        0.0, ek_xz, ek_h,
        0.0, ek_h, ek_yz;

    return Eigen::kroneckerProduct(ek, s0);
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
    rso << 0.0, 1i * std::sin(ky), 1i * std::sin(kx),
        -1i * std::sin(ky), 0.0, 0.0,
        -1i * std::sin(kx), 0.0, 0.0;
    Hamiltonian HRSO = delta_RSO * Eigen::kroneckerProduct(rso, s0);

    return HRSO;
}