#include "System2D.h"

#include <iostream>

void System2D::update_HBdG_nonzero_indices()
{
    HBdG_discrete_ky_onsite_nonzero_indices = assemble_HBdG_nonzero_indices(Hk_discrete_ky_onsite_nonzero_indices, Delta_discrete_ky_onsite_nonzero_indices, Delta_Adjoint_discrete_ky_onsite_nonzero_indices, mHmkT_discrete_ky_onsite_nonzero_indices);
    HBdG_discrete_ky_hopping_p_nonzero_indices = assemble_HBdG_nonzero_indices(Hk_discrete_ky_hopping_p_nonzero_indices, Delta_discrete_ky_hopping_p_nonzero_indices, Delta_Adjoint_discrete_ky_hopping_p_nonzero_indices, mHmkT_discrete_ky_hopping_p_nonzero_indices);
    HBdG_discrete_ky_hopping_m_nonzero_indices = assemble_HBdG_nonzero_indices(Hk_discrete_ky_hopping_m_nonzero_indices, Delta_discrete_ky_hopping_m_nonzero_indices, Delta_Adjoint_discrete_ky_hopping_m_nonzero_indices, mHmkT_discrete_ky_hopping_m_nonzero_indices);

    HBdG_discrete_onsite_nonzero_indices = assemble_HBdG_nonzero_indices(Hk_discrete_onsite_nonzero_indices, Delta_discrete_onsite_nonzero_indices, Delta_Adjoint_discrete_onsite_nonzero_indices, mHmkT_discrete_onsite_nonzero_indices);
    HBdG_discrete_hopping_xp_nonzero_indices = assemble_HBdG_nonzero_indices(Hk_discrete_hopping_xp_nonzero_indices, Delta_discrete_hopping_xp_nonzero_indices, Delta_Adjoint_discrete_hopping_xp_nonzero_indices, mHmkT_discrete_hopping_xp_nonzero_indices);
    HBdG_discrete_hopping_xm_nonzero_indices = assemble_HBdG_nonzero_indices(Hk_discrete_hopping_xm_nonzero_indices, Delta_discrete_hopping_xm_nonzero_indices, Delta_Adjoint_discrete_hopping_xm_nonzero_indices, mHmkT_discrete_hopping_xm_nonzero_indices);
    HBdG_discrete_hopping_yp_nonzero_indices = assemble_HBdG_nonzero_indices(Hk_discrete_hopping_yp_nonzero_indices, Delta_discrete_hopping_yp_nonzero_indices, Delta_Adjoint_discrete_hopping_yp_nonzero_indices, mHmkT_discrete_hopping_yp_nonzero_indices);
    HBdG_discrete_hopping_ym_nonzero_indices = assemble_HBdG_nonzero_indices(Hk_discrete_hopping_ym_nonzero_indices, Delta_discrete_hopping_ym_nonzero_indices, Delta_Adjoint_discrete_hopping_ym_nonzero_indices, mHmkT_discrete_hopping_ym_nonzero_indices);
}

Hamiltonian System2D::Hk_discrete_ky(double kx, std::size_t n_ky) const
{
    return assemble_matrix_discrete_ky(kx, n_ky, [this](double kx, double y)
                                       { return Hk_discrete_ky_onsite(kx, y); }, [this](double kx, double y)
                                       { return Hk_discrete_ky_hopping_p(kx, y); }, [this](double kx, double y)
                                       { return Hk_discrete_ky_hopping_m(kx, y); });
}

Hamiltonian System2D::Hk_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_matrix_discrete(n_kx, n_ky, [this](double x, double y)
                                    { return Hk_discrete_onsite(x, y); }, [this](double x, double y)
                                    { return Hk_discrete_hopping_xp(x, y); }, [this](double x, double y)
                                    { return Hk_discrete_hopping_xm(x, y); }, [this](double x, double y)
                                    { return Hk_discrete_hopping_yp(x, y); }, [this](double x, double y)
                                    { return Hk_discrete_hopping_ym(x, y); });
}

std::vector<Triplet> System2D::triplets_Hk_discrete_ky(double kx, std::size_t n_ky) const
{
    return assemble_triplets_discrete_ky(kx, n_ky, [this](double kx, double y)
                                         { return Hk_discrete_ky_onsite(kx, y); }, [this](double kx, double y)
                                         { return Hk_discrete_ky_hopping_p(kx, y); }, [this](double kx, double y)
                                         { return Hk_discrete_ky_hopping_m(kx, y); }, Hk_discrete_ky_onsite_nonzero_indices, Hk_discrete_ky_hopping_p_nonzero_indices, Hk_discrete_ky_hopping_m_nonzero_indices);
}

std::vector<Triplet> System2D::triplets_Hk_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_triplets_discrete(n_kx, n_ky, [this](double x, double y)
                                      { return Hk_discrete_onsite(x, y); }, [this](double x, double y)
                                      { return Hk_discrete_hopping_xp(x, y); }, [this](double x, double y)
                                      { return Hk_discrete_hopping_xm(x, y); }, [this](double x, double y)
                                      { return Hk_discrete_hopping_yp(x, y); }, [this](double x, double y)
                                      { return Hk_discrete_hopping_ym(x, y); }, Hk_discrete_onsite_nonzero_indices, Hk_discrete_hopping_xp_nonzero_indices, Hk_discrete_hopping_xm_nonzero_indices, Hk_discrete_hopping_yp_nonzero_indices, Hk_discrete_hopping_ym_nonzero_indices);
}

Hamiltonian System2D::HBdG(double kx, double ky) const
{
    return assemble_HBdG(Hk(kx, ky), Delta(kx, ky), Delta_Adjoint(kx, ky), mHmkT(kx, ky));
}

Hamiltonian System2D::HBdG_discrete_ky(double kx, std::size_t n_ky) const
{
    return assemble_matrix_discrete_ky(kx, n_ky, [this](double kx, double y)
                                       { return HBdG_discrete_ky_onsite(kx, y); }, [this](double kx, double y)
                                       { return HBdG_discrete_ky_hopping_p(kx, y); }, [this](double kx, double y)
                                       { return HBdG_discrete_ky_hopping_m(kx, y); });
}

Hamiltonian System2D::HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_matrix_discrete(n_kx, n_ky, [this](double x, double y)
                                    { return HBdG_discrete_onsite(x, y); }, [this](double x, double y)
                                    { return HBdG_discrete_hopping_xp(x, y); }, [this](double x, double y)
                                    { return HBdG_discrete_hopping_xm(x, y); }, [this](double x, double y)
                                    { return HBdG_discrete_hopping_yp(x, y); }, [this](double x, double y)
                                    { return HBdG_discrete_hopping_ym(x, y); });
}

std::vector<Triplet> System2D::triplets_HBdG_discrete_ky(double kx, std::size_t n_ky) const
{
    return assemble_triplets_discrete_ky(kx, n_ky, [this](double kx, double y)
                                         { return HBdG_discrete_ky_onsite(kx, y); }, [this](double kx, double y)
                                         { return HBdG_discrete_ky_hopping_p(kx, y); }, [this](double kx, double y)
                                         { return HBdG_discrete_ky_hopping_m(kx, y); }, HBdG_discrete_ky_onsite_nonzero_indices, HBdG_discrete_ky_hopping_p_nonzero_indices, HBdG_discrete_ky_hopping_m_nonzero_indices);
}

std::vector<Triplet> System2D::triplets_HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_triplets_discrete(n_kx, n_ky, [this](double x, double y)
                                      { return HBdG_discrete_onsite(x, y); }, [this](double x, double y)
                                      { return HBdG_discrete_hopping_xp(x, y); }, [this](double x, double y)
                                      { return HBdG_discrete_hopping_xm(x, y); }, [this](double x, double y)
                                      { return HBdG_discrete_hopping_yp(x, y); }, [this](double x, double y)
                                      { return HBdG_discrete_hopping_ym(x, y); }, HBdG_discrete_onsite_nonzero_indices, HBdG_discrete_hopping_xp_nonzero_indices, HBdG_discrete_hopping_xm_nonzero_indices, HBdG_discrete_hopping_yp_nonzero_indices, HBdG_discrete_hopping_ym_nonzero_indices);
}

Hamiltonian System2D::HBdG_discrete_ky_onsite(double kx, double y) const
{
    return assemble_HBdG(Hk_discrete_ky_onsite(kx, y), Delta_discrete_ky_onsite(kx, y), Delta_Adjoint_discrete_ky_onsite(kx, y), mHmkT_discrete_ky_onsite(kx, y));
}

Hamiltonian System2D::HBdG_discrete_ky_hopping_p(double kx, double y) const
{
    return assemble_HBdG(Hk_discrete_ky_hopping_p(kx, y), Delta_discrete_ky_hopping_p(kx, y), Delta_Adjoint_discrete_ky_hopping_p(kx, y), mHmkT_discrete_ky_hopping_p(kx, y));
}

Hamiltonian System2D::HBdG_discrete_ky_hopping_m(double kx, double y) const
{
    return assemble_HBdG(Hk_discrete_ky_hopping_m(kx, y), Delta_discrete_ky_hopping_m(kx, y), Delta_Adjoint_discrete_ky_hopping_m(kx, y), mHmkT_discrete_ky_hopping_m(kx, y));
}

Hamiltonian System2D::HBdG_discrete_onsite(double x, double y) const
{
    return assemble_HBdG(Hk_discrete_onsite(x, y), Delta_discrete_onsite(x, y), Delta_Adjoint_discrete_onsite(x, y), mHmkT_discrete_onsite(x, y));
}

Hamiltonian System2D::HBdG_discrete_hopping_xp(double x, double y) const
{
    return assemble_HBdG(Hk_discrete_hopping_xp(x, y), Delta_discrete_hopping_xp(x, y), Delta_Adjoint_discrete_hopping_xp(x, y), mHmkT_discrete_hopping_xp(x, y));
}

Hamiltonian System2D::HBdG_discrete_hopping_xm(double x, double y) const
{
    return assemble_HBdG(Hk_discrete_hopping_xm(x, y), Delta_discrete_hopping_xm(x, y), Delta_Adjoint_discrete_hopping_xm(x, y), mHmkT_discrete_hopping_xm(x, y));
}

Hamiltonian System2D::HBdG_discrete_hopping_yp(double x, double y) const
{
    return assemble_HBdG(Hk_discrete_hopping_yp(x, y), Delta_discrete_hopping_yp(x, y), Delta_Adjoint_discrete_hopping_yp(x, y), mHmkT_discrete_hopping_yp(x, y));
}

Hamiltonian System2D::HBdG_discrete_hopping_ym(double x, double y) const
{
    return assemble_HBdG(Hk_discrete_hopping_ym(x, y), Delta_discrete_hopping_ym(x, y), Delta_Adjoint_discrete_hopping_ym(x, y), mHmkT_discrete_hopping_ym(x, y));
}

Hamiltonian System2D::assemble_HBdG(const Hamiltonian &Hk_mat, const Hamiltonian &Delta_mat, const Hamiltonian &Delta_Adjoint_mat, const Hamiltonian &mHmkT_mat) const
{
    Hamiltonian H(2 * n_bands, 2 * n_bands);

    H << Hk_mat, Delta_mat,
        Delta_Adjoint_mat, mHmkT_mat;

    return H;
}

std::vector<ElementIndex> System2D::assemble_HBdG_nonzero_indices(const std::vector<ElementIndex> &Hk_nonzero_indices,
                                                                  const std::vector<ElementIndex> &Delta_nonzero_indices,
                                                                  const std::vector<ElementIndex> &Delta_Adjoint_nonzero_indices,
                                                                  const std::vector<ElementIndex> &mHmkT_nonzero_indices) const
{
    std::vector<ElementIndex> nonzero_indices;
    nonzero_indices.reserve(Hk_nonzero_indices.size() + Delta_nonzero_indices.size() + Delta_Adjoint_nonzero_indices.size() + mHmkT_nonzero_indices.size());

    for (const auto &t : Hk_nonzero_indices)
    {
        nonzero_indices.push_back(t);
    }

    for (const auto &t : Delta_nonzero_indices)
    {
        nonzero_indices.emplace_back(t.row, t.col + n_bands);
    }

    for (const auto &t : Delta_Adjoint_nonzero_indices)
    {
        nonzero_indices.emplace_back(t.row + n_bands, t.col);
    }

    for (const auto &t : mHmkT_nonzero_indices)
    {
        nonzero_indices.emplace_back(t.row + n_bands, t.col + n_bands);
    }

    return nonzero_indices;
}

std::vector<Triplet> System2D::assemble_triplets_HBdG(const std::vector<Triplet> &Hk_triplets, const std::vector<Triplet> &Delta_triplets, const std::vector<Triplet> &Delta_adjoint_triplets, const std::vector<Triplet> &mHmkT_triplets) const
{
    std::vector<Triplet> triplets;
    triplets.reserve(Hk_triplets.size() + Delta_triplets.size() + Delta_adjoint_triplets.size() + mHmkT_triplets.size());

    for (const auto &t : Hk_triplets)
    {
        triplets.push_back(t);
    }

    for (const auto &t : Delta_triplets)
    {
        triplets.emplace_back(t.row(), t.col() + n_bands, t.value());
    }

    for (const auto &t : Delta_adjoint_triplets)
    {
        triplets.emplace_back(t.row() + n_bands, t.col(), t.value());
    }

    for (const auto &t : mHmkT_triplets)
    {
        triplets.emplace_back(t.row() + n_bands, t.col() + n_bands, t.value());
    }

    return triplets;
}

Hamiltonian System2D::assemble_matrix_discrete_ky(double kx, std::size_t n_ky,
                                                  const HamiltonianFunction &onsite,
                                                  const HamiltonianFunction &hopping_p,
                                                  const HamiltonianFunction &hopping_m) const
{
    std::size_t submatrix_size = onsite(kx, 0).rows();
    Hamiltonian H = Hamiltonian::Zero(submatrix_size * n_ky, submatrix_size * n_ky);
    int iy = 0;

    for (int i = 0; i < n_ky; ++i)
    {
        double y = dy * i;

        Hamiltonian o = onsite(kx, y);
        Hamiltonian hp = hopping_p(kx, y);
        Hamiltonian hm = hopping_m(kx, y);

        H.block(i * submatrix_size, i * submatrix_size, submatrix_size, submatrix_size) = o;

        iy = i - 1;
        if (i % n_ky != 0)
            H.block(i * submatrix_size, iy * submatrix_size, submatrix_size, submatrix_size) = hm;

        iy = i + 1;
        if (iy % n_ky != 0)
            H.block(i * submatrix_size, iy * submatrix_size, submatrix_size, submatrix_size) = hp;
    }

    return H;
}

Hamiltonian System2D::assemble_matrix_discrete(std::size_t n_kx, std::size_t n_ky,
                                               const HamiltonianFunction &onsite,
                                               const HamiltonianFunction &hopping_xp,
                                               const HamiltonianFunction &hopping_xm,
                                               const HamiltonianFunction &hopping_yp,
                                               const HamiltonianFunction &hopping_ym) const
{
    std::size_t submatrix_size = onsite(0, 0).rows();
    Hamiltonian H = Hamiltonian::Zero(submatrix_size * n_kx * n_ky, submatrix_size * n_kx * n_ky);
    int ix = 0;
    int iy = 0;

    for (int i = 0; i < n_kx * n_ky; ++i)
    {
        double x = dx * (i / n_ky);
        double y = dy * (i % n_ky);

        Hamiltonian o = onsite(x, y);
        Hamiltonian hxp = hopping_xp(x, y);
        Hamiltonian hxm = hopping_xm(x, y);
        Hamiltonian hyp = hopping_yp(x, y);
        Hamiltonian hym = hopping_ym(x, y);

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
    }

    return H;
}

std::vector<Triplet> System2D::assemble_triplets_discrete_ky(double kx, std::size_t n_ky,
                                                             const HamiltonianFunction &onsite,
                                                             const HamiltonianFunction &hopping_p,
                                                             const HamiltonianFunction &hopping_m,
                                                             const std::vector<ElementIndex> &onsite_nonzero_indices,
                                                             const std::vector<ElementIndex> &hopping_p_nonzero_indices,
                                                             const std::vector<ElementIndex> &hopping_m_nonzero_indices) const
{
    std::size_t submatrix_size = onsite(kx, 0).rows();
    std::size_t n_elements = onsite_nonzero_indices.size() + hopping_p_nonzero_indices.size() + hopping_m_nonzero_indices.size();
    std::vector<Triplet> triplets;
    triplets.reserve(n_elements * n_ky);
    int iy = 0;

    for (int i = 0; i < n_ky; ++i)
    {
        double y = dy * i;

        for (const auto &index : onsite_nonzero_indices)
        {
            auto [row, col] = index;
            triplets.emplace_back(i * submatrix_size + row, i * submatrix_size + col, onsite(kx, y)(row, col));
        }

        iy = i - 1;
        if (i % n_ky != 0)
            for (const auto &index : hopping_m_nonzero_indices)
            {
                auto [row, col] = index;
                triplets.emplace_back(i * submatrix_size + row, iy * submatrix_size + col, hopping_m(kx, y)(row, col));
            }

        iy = i + 1;
        if (iy % n_ky != 0)
            for (const auto &index : hopping_p_nonzero_indices)
            {
                auto [row, col] = index;
                triplets.emplace_back(i * submatrix_size + row, iy * submatrix_size + col, hopping_p(kx, y)(row, col));
            }
    }

    return triplets;
}

std::vector<Triplet> System2D::assemble_triplets_discrete(std::size_t n_kx, std::size_t n_ky,
                                                          const HamiltonianFunction &onsite,
                                                          const HamiltonianFunction &hopping_xp,
                                                          const HamiltonianFunction &hopping_xm,
                                                          const HamiltonianFunction &hopping_yp,
                                                          const HamiltonianFunction &hopping_ym,
                                                          const std::vector<ElementIndex> &onsite_nonzero_indices,
                                                          const std::vector<ElementIndex> &hopping_xp_nonzero_indices,
                                                          const std::vector<ElementIndex> &hopping_xm_nonzero_indices,
                                                          const std::vector<ElementIndex> &hopping_yp_nonzero_indices,
                                                          const std::vector<ElementIndex> &hopping_ym_nonzero_indices) const
{
    std::size_t submatrix_size = onsite(0, 0).rows();
    std::size_t n_elements = onsite_nonzero_indices.size() + hopping_xp_nonzero_indices.size() + hopping_xm_nonzero_indices.size() +
                             hopping_yp_nonzero_indices.size() + hopping_ym_nonzero_indices.size();
    std::vector<Triplet> triplets;
    triplets.reserve(n_elements * n_kx * n_ky);
    int ix = 0;
    int iy = 0;

    for (int i = 0; i < n_kx * n_ky; ++i)
    {
        double x = dx * (i / n_ky);
        double y = dy * (i % n_ky);

        for (const auto &index : onsite_nonzero_indices)
        {
            auto [row, col] = index;
            triplets.emplace_back(i * submatrix_size + row, i * submatrix_size + col, onsite(x, y)(row, col));
        }

        iy = i - 1;
        if (i % n_ky != 0)
            for (const auto &index : hopping_ym_nonzero_indices)
            {
                auto [row, col] = index;
                triplets.emplace_back(i * submatrix_size + row, iy * submatrix_size + col, hopping_ym(x, y)(row, col));
            }

        iy = i + 1;
        if (iy % n_ky != 0)
            for (const auto &index : hopping_yp_nonzero_indices)
            {
                auto [row, col] = index;
                triplets.emplace_back(i * submatrix_size + row, iy * submatrix_size + col, hopping_yp(x, y)(row, col));
            }

        ix = i - n_ky;
        if (ix >= 0)
            for (const auto &index : hopping_xm_nonzero_indices)
            {
                auto [row, col] = index;
                triplets.emplace_back(i * submatrix_size + row, ix * submatrix_size + col, hopping_xm(x, y)(row, col));
            }

        ix = i + n_ky;
        if (ix < n_kx * n_ky)
            for (const auto &index : hopping_xp_nonzero_indices)
            {
                auto [row, col] = index;
                triplets.emplace_back(i * submatrix_size + row, ix * submatrix_size + col, hopping_xp(x, y)(row, col));
            }
    }

    return triplets;
}
