#include "System2D.h"

Hamiltonian System2D::Hk(double kx, double ky) const
{
    return assemble_matrix(Hk_triplets(kx, ky), n_bands, n_bands);
}

Hamiltonian System2D::Hk_discrete_ky(double kx, std::size_t n_ky) const
{
    return assemble_matrix(Hk_discrete_ky_triplets(kx, n_ky), n_bands * n_ky, n_bands * n_ky);
}

SparseHamiltonian System2D::Hk_discrete_ky_sparse(double kx, std::size_t n_ky) const
{
    return assemble_sparse_matrix(Hk_discrete_ky_triplets(kx, n_ky), n_bands * n_ky, n_bands * n_ky);
}

Hamiltonian System2D::Hk_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_matrix(Hk_discrete_triplets(n_kx, n_ky), n_bands * n_kx * n_ky, n_bands * n_kx * n_ky);
}

SparseHamiltonian System2D::Hk_discrete_sparse(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_sparse_matrix(Hk_discrete_triplets(n_kx, n_ky), n_bands * n_kx * n_ky, n_bands * n_kx * n_ky);
}

Hamiltonian System2D::HBdG(double kx, double ky) const
{
    return assemble_matrix(HBdG_triplets(kx, ky), n_bands * 2, n_bands * 2);
}

Hamiltonian System2D::HBdG_discrete_ky(double kx, std::size_t n_ky) const
{
    return assemble_matrix(HBdG_discrete_ky_triplets(kx, n_ky), n_bands * 2 * n_ky, n_bands * 2 * n_ky);
}

SparseHamiltonian System2D::HBdG_discrete_ky_sparse(double kx, std::size_t n_ky) const
{
    return assemble_sparse_matrix(HBdG_discrete_ky_triplets(kx, n_ky), n_bands * 2 * n_ky, n_bands * 2 * n_ky);
}

Hamiltonian System2D::HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_matrix(HBdG_discrete_triplets(n_kx, n_ky), n_bands * 2 * n_kx * n_ky, n_bands * 2 * n_kx * n_ky);
}

SparseHamiltonian System2D::HBdG_discrete_sparse(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_sparse_matrix(HBdG_discrete_triplets(n_kx, n_ky), n_bands * 2 * n_kx * n_ky, n_bands * 2 * n_kx * n_ky);
}

std::vector<Triplet> System2D::generate_triplets(const Hamiltonian &H, bool only_upper_triangular) const
{
    std::vector<Triplet> triplets;
    for (std::size_t i = 0; i < H.rows(); ++i)
        for (std::size_t j = only_upper_triangular ? i : 0; j < H.cols(); ++j)
            triplets.emplace_back(i, j, H(i, j));
    return triplets;
}

std::vector<Triplet> System2D::negate_triplets(const std::vector<Triplet> &triplets) const
{
    std::vector<Triplet> negated_triplets;
    negated_triplets.reserve(triplets.size());
    for (const auto &tr : triplets)
    {
        negated_triplets.emplace_back(tr.row(), tr.col(), -tr.value());
    }
    return negated_triplets;
}

std::vector<Triplet> System2D::Hk_discrete_ky_triplets(double kx, std::size_t n_ky) const
{
    std::size_t submatrix_size = n_bands;

    std::vector<Triplet> triplets = assemble_triplets_discrete_ky(
        [this](double kx, double y)
        { return Hk_discrete_ky_onsite_triplets(kx, y); },
        [this](double kx, double y)
        { return Hk_discrete_ky_hopping_p_triplets(kx, y); },
        kx, n_ky, submatrix_size);

    return triplets;
}

std::vector<Triplet> System2D::Hk_discrete_triplets(std::size_t n_kx, std::size_t n_ky) const
{
    std::size_t submatrix_size = n_bands;

    std::vector<Triplet> triplets = assemble_triplets_discrete(
        [this](double x, double y)
        { return Hk_discrete_onsite_triplets(x, y); },
        [this](double x, double y)
        { return Hk_discrete_hopping_xp_triplets(x, y); },
        [this](double x, double y)
        { return Hk_discrete_hopping_yp_triplets(x, y); },
        [this](double x, double y)
        { return Hk_discrete_hopping_pp_triplets(x, y); },
        [this](double x, double y)
        { return Hk_discrete_hopping_pm_triplets(x, y); },
        n_kx, n_ky, submatrix_size);

    return triplets;
}

std::vector<Triplet> System2D::join_triplets_for_HBdG(
    const std::vector<Triplet> &Hk_tr,
    const std::vector<Triplet> &Delta_tr,
    const std::vector<Triplet> &mHmkT_tr) const
{
    std::vector<Triplet> triplets = Hk_tr;
    triplets.reserve(Hk_tr.size() + Delta_tr.size() + mHmkT_tr.size());

    for (const auto &tr : Delta_tr)
    {
        triplets.emplace_back(tr.row(), tr.col() + n_bands, tr.value());
    }

    for (const auto &tr : mHmkT_tr)
    {
        triplets.emplace_back(tr.row() + n_bands, tr.col() + n_bands, tr.value());
    }

    return triplets;
}

std::vector<Triplet> System2D::HBdG_triplets(double kx, double ky) const
{
    std::vector<Triplet> Hk_tr = Hk_triplets(kx, ky);
    std::vector<Triplet> Delta_tr = Delta_triplets(kx, ky);
    std::vector<Triplet> mHmkT_tr = mHmkT_triplets(kx, ky);

    return join_triplets_for_HBdG(Hk_tr, Delta_tr, mHmkT_tr);
}

std::vector<Triplet> System2D::HBdG_discrete_ky_triplets(double kx, std::size_t n_ky) const
{
    std::size_t submatrix_size = n_bands * 2;

    std::vector<Triplet> Hk_tr = assemble_triplets_discrete_ky(
        [this](double kx, double y)
        { return Hk_discrete_ky_onsite_triplets(kx, y); },
        [this](double kx, double y)
        { return Hk_discrete_ky_hopping_p_triplets(kx, y); },
        kx, n_ky, submatrix_size);

    std::vector<Triplet> Delta_tr = assemble_triplets_discrete_ky(
        [this](double kx, double y)
        { return Delta_discrete_ky_onsite_triplets(kx, y); },
        [this](double kx, double y)
        { return Delta_discrete_ky_hopping_p_triplets(kx, y); },
        kx, n_ky, submatrix_size);

    std::vector<Triplet> mHmkT_tr = assemble_triplets_discrete_ky(
        [this](double kx, double y)
        { return mHmkT_discrete_ky_onsite_triplets(kx, y); },
        [this](double kx, double y)
        { return mHmkT_discrete_ky_hopping_p_triplets(kx, y); },
        kx, n_ky, submatrix_size);

    return join_triplets_for_HBdG(Hk_tr, Delta_tr, mHmkT_tr);
}

std::vector<Triplet> System2D::HBdG_discrete_triplets(std::size_t n_kx, std::size_t n_ky) const
{
    std::size_t submatrix_size = n_bands * 2;

    std::vector<Triplet> Hk_tr = assemble_triplets_discrete(
        [this](double x, double y)
        { return Hk_discrete_onsite_triplets(x, y); },
        [this](double x, double y)
        { return Hk_discrete_hopping_xp_triplets(x, y); },
        [this](double x, double y)
        { return Hk_discrete_hopping_yp_triplets(x, y); },
        [this](double x, double y)
        { return Hk_discrete_hopping_pp_triplets(x, y); },
        [this](double x, double y)
        { return Hk_discrete_hopping_pm_triplets(x, y); },
        n_kx, n_ky, submatrix_size);

    std::vector<Triplet> Delta_tr = assemble_triplets_discrete(
        [this](double x, double y)
        { return Delta_discrete_onsite_triplets(x, y); },
        [this](double x, double y)
        { return Delta_discrete_hopping_xp_triplets(x, y); },
        [this](double x, double y)
        { return Delta_discrete_hopping_yp_triplets(x, y); },
        [this](double x, double y)
        { return Delta_discrete_hopping_pp_triplets(x, y); },
        [this](double x, double y)
        { return Delta_discrete_hopping_pm_triplets(x, y); },
        n_kx, n_ky, submatrix_size);
    std::vector<Triplet> mHmkT_tr = assemble_triplets_discrete(
        [this](double x, double y)
        { return mHmkT_discrete_onsite_triplets(x, y); },
        [this](double x, double y)
        { return mHmkT_discrete_hopping_xp_triplets(x, y); },
        [this](double x, double y)
        { return mHmkT_discrete_hopping_yp_triplets(x, y); },
        [this](double x, double y)
        { return mHmkT_discrete_hopping_pp_triplets(x, y); },
        [this](double x, double y)
        { return mHmkT_discrete_hopping_pm_triplets(x, y); },
        n_kx, n_ky, submatrix_size);

    return join_triplets_for_HBdG(Hk_tr, Delta_tr, mHmkT_tr);
}

void System2D::append_triplets(std::size_t row_offset, std::size_t col_offset,
                               std::vector<Triplet> &target,
                               const std::vector<Triplet> &source) const
{
    for (const auto &t : source)
    {
        target.emplace_back(t.row() + row_offset, t.col() + col_offset, t.value());
    }
}

std::vector<Triplet> System2D::assemble_triplets_discrete_ky(const TripletFunc &onsite_tf,
                                                             const TripletFunc &hopping_p_tf,
                                                             double kx, std::size_t n_ky, std::size_t submatrix_size) const
{
    std::vector<Triplet> triplets;
    triplets.reserve(2 * n_ky * (onsite_tf(kx, 0.0).size() + hopping_p_tf(kx, 0.0).size()));

    double y = 0.0;

    int is = 0;

    int j = 0;
    int js = 0;

    for (int i = 0; i < n_ky; ++i)
    {
        y = dy * i;

        is = i * submatrix_size;

        // osite
        j = i;
        js = j * submatrix_size;
        append_triplets(is, js, triplets, onsite_tf(kx, y));

        // hopping_p
        j = i + 1;
        js = j * submatrix_size;
        if (j % n_ky != 0)
            append_triplets(is, js, triplets, hopping_p_tf(kx, y));
    }

    return triplets;
}

std::vector<Triplet> System2D::assemble_triplets_discrete(const TripletFunc &onsite_tf,
                                                          const TripletFunc &hopping_xp_tf,
                                                          const TripletFunc &hopping_yp_tf,
                                                          const TripletFunc &hopping_pp_tf,
                                                          const TripletFunc &hopping_pm_tf,
                                                          std::size_t n_kx, std::size_t n_ky, std::size_t submatrix_size) const
{
    std::vector<Triplet> triplets;
    triplets.reserve(2 * n_kx * n_ky * (onsite_tf(0.0, 0.0).size() + hopping_xp_tf(0.0, 0.0).size() + hopping_yp_tf(0.0, 0.0).size() + hopping_pp_tf(0.0, 0.0).size() + hopping_pm_tf(0.0, 0.0).size()));

    double x = 0.0;
    double y = 0.0;

    int is = 0;

    int j = 0;
    int js = 0;

    for (int i = 0; i < n_kx; ++i)
    {
        x = dx * i;
        append_triplets(i * n_ky * submatrix_size, i * n_ky * submatrix_size, triplets, assemble_triplets_discrete_ky(onsite_tf, hopping_yp_tf, x, n_ky, submatrix_size));
    }

    for (int i = 0; i < n_kx * n_ky; ++i)
    {
        x = dx * (i / n_ky);
        y = dy * (i % n_ky);

        is = i * submatrix_size;

        // hopping_xp
        j = i + n_ky;
        js = j * submatrix_size;
        if (j < n_kx * n_ky)
            append_triplets(is, js, triplets, hopping_xp_tf(x, y));

        // hopping_pp
        j = i + n_ky + 1;
        js = j * submatrix_size;
        if (j % n_ky != 0 && j < n_kx * n_ky)
            append_triplets(is, js, triplets, hopping_pp_tf(x, y));

        // hopping_pm
        j = i + n_ky - 1;
        js = j * submatrix_size;
        if ((j + 1) % n_ky != 0 && j < n_kx * n_ky)
            append_triplets(is, js, triplets, hopping_pm_tf(x, y));
    }

    return triplets;
}

Hamiltonian System2D::assemble_matrix(const std::vector<Triplet> &triplets, std::size_t n_rows, std::size_t n_cols) const
{
    Hamiltonian H = Hamiltonian::Zero(n_rows, n_cols);

    for (const auto &t : triplets)
    {
        H(t.row(), t.col()) = t.value();

        if (t.row() != t.col()) // add conjugate for off-diagonal elements
        {
            H(t.col(), t.row()) = std::conj(t.value());
        }
    }

    return H;
}

SparseHamiltonian System2D::assemble_sparse_matrix(const std::vector<Triplet> &triplets, std::size_t n_rows, std::size_t n_cols) const
{
    SparseHamiltonian H(n_rows, n_cols);

    std::vector<Triplet> full_triplets = triplets;
    full_triplets.reserve(triplets.size() * 2);

    for (const auto &t : triplets)
    {
        if (t.row() != t.col()) // add conjugate for off-diagonal elements
        {
            full_triplets.emplace_back(t.col(), t.row(), std::conj(t.value()));
        }
    }

    H.setFromTriplets(full_triplets.begin(), full_triplets.end());

    return H;
}