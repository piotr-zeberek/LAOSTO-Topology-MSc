#include "System2DMatrixBased.h"
#include "utils.h"

#include <vector>

Hamiltonian System2DMatrixBased::HBdG_discrete_ky(double kx, std::size_t n_ky) const
{
    Hamiltonian H = Hamiltonian::Zero(2 * n_bands * n_ky, 2 * n_bands * n_ky);

    int iy = 0;

    for (int i = 0; i < n_ky; ++i)
    {
        double y = dy * i;

        Hamiltonian o = HBdG_discrete_ky_onsite(kx, y);
        Hamiltonian hp = HBdG_discrete_ky_hopping_p(kx, y);
        Hamiltonian hm = HBdG_discrete_ky_hopping_m(kx, y);

        H.block(i * 2 * n_bands, i * 2 * n_bands, 2 * n_bands, 2 * n_bands) = o;

        iy = i - 1;
        if (i % n_ky != 0)
            H.block(i * 2 * n_bands, iy * 2 * n_bands, 2 * n_bands, 2 * n_bands) = hm;

        iy = i + 1;
        if (iy % n_ky != 0)
            H.block(i * 2 * n_bands, iy * 2 * n_bands, 2 * n_bands, 2 * n_bands) = hp;
    }

    return H;
}

Triplets System2DMatrixBased::triplets_HBdG_discrete_ky(double kx, std::size_t n_ky) const
{

    Triplets triplets;
    triplets.reserve(2 * n_bands * 2 * n_bands * n_ky);

    int iy = 0;

    for (int i = 0; i < n_ky; ++i)
    {
        double y = dy * i;

        SparseHamiltonian o = HBdG_discrete_ky_onsite(kx, y).sparseView();
        SparseHamiltonian hp = HBdG_discrete_ky_hopping_p(kx, y).sparseView();
        SparseHamiltonian hm = HBdG_discrete_ky_hopping_m(kx, y).sparseView();

        Triplets o_triplets = get_triplets(o);
        Triplets hp_triplets = get_triplets(hp);
        Triplets hm_triplets = get_triplets(hm);

        add_triplets(triplets, o_triplets, i * 2 * n_bands, i * 2 * n_bands);

        iy = i - 1;
        if (i % n_ky != 0)
        {
            add_triplets(triplets, hm_triplets, i * 2 * n_bands, iy * 2 * n_bands);
        }

        iy = i + 1;
        if (iy % n_ky != 0)
        {
            add_triplets(triplets, hp_triplets, i * 2 * n_bands, iy * 2 * n_bands);
        }
    }

    return triplets;
}

Hamiltonian System2DMatrixBased::HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    Hamiltonian H = Hamiltonian::Zero(2 * n_bands * n_kx * n_ky, 2 * n_bands * n_kx * n_ky);

    int ix = 0;
    int iy = 0;

    for (int i = 0; i < n_kx * n_ky; ++i)
    {
        double x = dx * (i % n_kx);
        double y = dy * (i / n_kx);

        Hamiltonian o = HBdG_discrete_onsite(x, y);
        Hamiltonian hxp = HBdG_discrete_hopping_xp(x, y);
        Hamiltonian hxm = HBdG_discrete_hopping_xm(x, y);
        Hamiltonian hyp = HBdG_discrete_hopping_yp(x, y);
        Hamiltonian hym = HBdG_discrete_hopping_ym(x, y);

        H.block(i * 2 * n_bands, i * 2 * n_bands, 2 * n_bands, 2 * n_bands) = o;

        iy = i - 1;
        if (i % n_ky != 0)
            H.block(i * 2 * n_bands, iy * 2 * n_bands, 2 * n_bands, 2 * n_bands) = hym;

        iy = i + 1;
        if (iy % n_ky != 0)
            H.block(i * 2 * n_bands, iy * 2 * n_bands, 2 * n_bands, 2 * n_bands) = hyp;

        ix = i - n_ky;
        if (ix >= 0)
            H.block(i * 2 * n_bands, ix * 2 * n_bands, 2 * n_bands, 2 * n_bands) = hxm;

        ix = i + n_ky;
        if (ix < n_kx * n_ky)
            H.block(i * 2 * n_bands, ix * 2 * n_bands, 2 * n_bands, 2 * n_bands) = hxp;
    }

    return H;
}

Triplets System2DMatrixBased::triplets_HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    Triplets triplets;
    triplets.reserve(2 * n_bands * 2 * n_bands * n_kx * n_ky);

    int ix = 0;
    int iy = 0;

    for (int i = 0; i < n_kx * n_ky; ++i)
    {
        double x = dx * (i % n_kx);
        double y = dy * (i / n_kx);

        SparseHamiltonian o = HBdG_discrete_onsite(x, y).sparseView();
        SparseHamiltonian hxp = HBdG_discrete_hopping_xp(x, y).sparseView();
        SparseHamiltonian hxm = HBdG_discrete_hopping_xm(x, y).sparseView();
        SparseHamiltonian hyp = HBdG_discrete_hopping_yp(x, y).sparseView();
        SparseHamiltonian hym = HBdG_discrete_hopping_ym(x, y).sparseView();

        Triplets o_triplets = get_triplets(o);
        Triplets hxp_triplets = get_triplets(hxp);
        Triplets hxm_triplets = get_triplets(hxm);
        Triplets hyp_triplets = get_triplets(hyp);
        Triplets hym_triplets = get_triplets(hym);

        add_triplets(triplets, o_triplets, i * 2 * n_bands, i * 2 * n_bands);

        iy = i - 1;
        if (i % n_ky != 0)
        {
            add_triplets(triplets, hym_triplets, i * 2 * n_bands, iy * 2 * n_bands);
        }

        iy = i + 1;
        if (iy % n_ky != 0)
        {
            add_triplets(triplets, hyp_triplets, i * 2 * n_bands, iy * 2 * n_bands);
        }

        ix = i - n_ky;
        if (ix >= 0)
        {
            add_triplets(triplets, hxm_triplets, i * 2 * n_bands, ix * 2 * n_bands);
        }

        ix = i + n_ky;
        if (ix < n_kx * n_ky)
        {
            add_triplets(triplets, hxp_triplets, i * 2 * n_bands, ix * 2 * n_bands);
        }
    }

    return triplets;
}

Triplets System2DMatrixBased::get_triplets(const SparseHamiltonian &H, double tol) const
{
    Triplets triplets;
    triplets.reserve(H.nonZeros());

    for (int k = 0; k < H.outerSize(); ++k)
    {
        for (SparseHamiltonian::InnerIterator it(H, k); it; ++it)
        {
            if (std::abs(it.value()) > tol)
            {
                triplets.push_back(Triplet(it.row(), it.col(), it.value()));
            }
        }
    }

    return triplets;
}

void System2DMatrixBased::add_triplets(Triplets &triplets, const Triplets &triplets_to_add, int row_offset, int col_offset) const
{
    for (const auto &t : triplets_to_add)
    {
        triplets.push_back(Triplet(t.row() + row_offset, t.col() + col_offset, t.value()));
    }
}