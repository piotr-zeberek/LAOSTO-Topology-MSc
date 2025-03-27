#include "ToyModel.h"

#include <unsupported/Eigen/KroneckerProduct>

Hamiltonian ToyModel::Hk(double kx, double ky) const
{
    return Hkin(kx, ky) + HZeeman(kx, ky) + HRashba(kx, ky) - mu * s0;
}

Hamiltonian ToyModel::HBdG(double kx, double ky) const
{
    Hamiltonian res = -delta_SC * Eigen::kroneckerProduct(sy, sy);

    res.topLeftCorner(2, 2) = Hk(kx, ky);
    res.bottomRightCorner(2, 2) = -Hk(-kx, -ky).transpose();

    return res;
}

Hamiltonian ToyModel::HBdG_discrete_ky_onsite(double kx, double y) const
{
    Hamiltonian res = -delta_SC * Eigen::kroneckerProduct(sy, sy);

    res.topLeftCorner(2, 2) = (2.0 * tx * (1.0 - std::cos(kx)) + 2.0 * ty - mu) * s0 + HZeeman(kx, 0.0) - delta_RSO_x * sy * std::sin(kx);
    res.bottomRightCorner(2, 2) = -(2.0 * tx * (1.0 - std::cos(kx)) + 2.0 * ty - mu) * s0 - HZeeman(kx, 0.0).transpose() + delta_RSO_x * sy * std::sin(kx);

    return res;
}

Hamiltonian ToyModel::HBdG_discrete_ky_hopping_p(double kx, double y) const
{
    Hamiltonian res = Hamiltonian::Zero(4, 4);

    res.topLeftCorner(2, 2) = -ty * s0 - 1i / 2.0 * delta_RSO_y * sx;
    res.bottomRightCorner(2, 2) = ty * s0 - 1i / 2.0 * delta_RSO_y * sx;

    return res;
}

Hamiltonian ToyModel::HBdG_discrete_ky_hopping_m(double kx, double y) const
{
    Hamiltonian res = Hamiltonian::Zero(4, 4);

    res.topLeftCorner(2, 2) = -ty * s0 + 1i / 2.0 * delta_RSO_y * sx;
    res.bottomRightCorner(2, 2) = ty * s0 + 1i / 2.0 * delta_RSO_y * sx;

    return res;
}

Hamiltonian ToyModel::HBdG_discrete_onsite(double x, double y) const
{
    Hamiltonian res = -delta_SC * Eigen::kroneckerProduct(sy, sy);

    double ek = 2.0 * (tx + ty);

    res.topLeftCorner(2, 2) = (ek - mu) * s0 + HZeeman(0.0, 0.0);
    res.bottomRightCorner(2, 2) = -res.topLeftCorner(2, 2).transpose();

    return res;
}

Hamiltonian ToyModel::HBdG_discrete_hopping_xp(double x, double y) const
{
    Hamiltonian res = Hamiltonian::Zero(4, 4);

    res.topLeftCorner(2, 2) = -tx * s0 + 1i * delta_RSO_x / 2.0 * sy;
    res.bottomRightCorner(2, 2) = tx * s0 - 1i * delta_RSO_x / 2.0 * sy;

    return res;
}

Hamiltonian ToyModel::HBdG_discrete_hopping_xm(double x, double y) const
{
    Hamiltonian res = Hamiltonian::Zero(4, 4);

    res.topLeftCorner(2, 2) = -tx * s0 - 1i * delta_RSO_x / 2.0 * sy;
    res.bottomRightCorner(2, 2) = tx * s0 + 1i * delta_RSO_x / 2.0 * sy;

    return res;
}

Hamiltonian ToyModel::HBdG_discrete_hopping_yp(double x, double y) const
{
    Hamiltonian res = Hamiltonian::Zero(4, 4);

    res.topLeftCorner(2, 2) = -ty * s0 - 1i * delta_RSO_y / 2.0 * sx;
    res.bottomRightCorner(2, 2) = ty * s0 - 1i * delta_RSO_y / 2.0 * sx;

    return res;
}

Hamiltonian ToyModel::HBdG_discrete_hopping_ym(double x, double y) const
{
    Hamiltonian res = Hamiltonian::Zero(4, 4);

    res.topLeftCorner(2, 2) = -ty * s0 + 1i * delta_RSO_y / 2.0 * sx;
    res.bottomRightCorner(2, 2) = ty * s0 + 1i * delta_RSO_y / 2.0 * sx;

    return res;
}

Hamiltonian ToyModel::Hkin(double kx, double ky) const
{
    return 2.0 * (tx * (1.0 - std::cos(kx)) + ty * (1.0 - std::cos(ky))) * s0;
}
Hamiltonian ToyModel::HZeeman(double kx, double ky) const
{
    return 0.5 * g_Lande * 0.5 * (Bx * sx + By * sy + Bz * sz);
}
Hamiltonian ToyModel::HRashba(double kx, double ky) const
{
    return delta_RSO_y * std::sin(ky) * sx - delta_RSO_x * std::sin(kx) * sy;
}