#include "System2D.h"
#include "utils.h"

#include <vector>

Hamiltonian System2D::HBdG_discrete_ky(double kx, std::size_t n_ky) const
{
    Hamiltonian H = Hamiltonian::Zero(n_bands_sc * n_ky, n_bands_sc * n_ky);

    int iy = 0;

    for (int i = 0; i < n_ky; ++i)
    {
        double y = dy * i;

        Hamiltonian o = HBdG_discrete_ky_onsite(kx, y);
        Hamiltonian hp = HBdG_discrete_ky_hopping_p(kx, y);
        Hamiltonian hm = HBdG_discrete_ky_hopping_m(kx, y);

        H.block(i * n_bands_sc, i * n_bands_sc, n_bands_sc, n_bands_sc) = o;

        iy = i - 1;
        if (i % n_ky != 0)
            H.block(i * n_bands_sc, iy * n_bands_sc, n_bands_sc, n_bands_sc) = hm;

        iy = i + 1;
        if (iy % n_ky != 0)
            H.block(i * n_bands_sc, iy * n_bands_sc, n_bands_sc, n_bands_sc) = hp;
    }

    return H;
}

std::vector<Triplet> System2D::HBdG_discrete_ky_triplets(double kx, std::size_t n_ky) const
{

    std::vector<Triplet> triplets;
    triplets.reserve(n_bands_sc * n_bands_sc * n_ky);

    int iy = 0;

    for (int i = 0; i < n_ky; ++i)
    {
        double y = dy * i;

        SparseHamiltonian o = HBdG_discrete_ky_onsite(kx, y).sparseView();
        SparseHamiltonian hp = HBdG_discrete_ky_hopping_p(kx, y).sparseView();
        SparseHamiltonian hm = HBdG_discrete_ky_hopping_m(kx, y).sparseView();

        std::vector<Triplet> o_triplets = get_triplets(o);
        std::vector<Triplet> hp_triplets = get_triplets(hp);
        std::vector<Triplet> hm_triplets = get_triplets(hm);

        add_triplets(triplets, o_triplets, i * n_bands_sc, i * n_bands_sc);

        iy = i - 1;
        if (i % n_ky != 0)
        {
            add_triplets(triplets, hm_triplets, i * n_bands_sc, iy * n_bands_sc);
        }

        iy = i + 1;
        if (iy % n_ky != 0)
        {
            add_triplets(triplets, hp_triplets, i * n_bands_sc, iy * n_bands_sc);
        }
    }

    return triplets;
}

Hamiltonian System2D::HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    Hamiltonian H = Hamiltonian::Zero(n_bands_sc * n_kx * n_ky, n_bands_sc * n_kx * n_ky);

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

        H.block(i * n_bands_sc, i * n_bands_sc, n_bands_sc, n_bands_sc) = o;

        iy = i - 1;
        if (i % n_ky != 0)
            H.block(i * n_bands_sc, iy * n_bands_sc, n_bands_sc, n_bands_sc) = hym;

        iy = i + 1;
        if (iy % n_ky != 0)
            H.block(i * n_bands_sc, iy * n_bands_sc, n_bands_sc, n_bands_sc) = hyp;

        ix = i - n_ky;
        if (ix >= 0)
            H.block(i * n_bands_sc, ix * n_bands_sc, n_bands_sc, n_bands_sc) = hxm;

        ix = i + n_ky;
        if (ix < n_kx * n_ky)
            H.block(i * n_bands_sc, ix * n_bands_sc, n_bands_sc, n_bands_sc) = hxp;
    }

    return H;
}

std::vector<Triplet> System2D::HBdG_discrete_triplets(std::size_t n_kx, std::size_t n_ky) const
{
    std::vector<Triplet> triplets;
    triplets.reserve(n_bands_sc * n_bands_sc * n_kx * n_ky);

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

        std::vector<Triplet> o_triplets = get_triplets(o);
        std::vector<Triplet> hxp_triplets = get_triplets(hxp);
        std::vector<Triplet> hxm_triplets = get_triplets(hxm);
        std::vector<Triplet> hyp_triplets = get_triplets(hyp);
        std::vector<Triplet> hym_triplets = get_triplets(hym);

        add_triplets(triplets, o_triplets, i * n_bands_sc, i * n_bands_sc);

        iy = i - 1;
        if (i % n_ky != 0)
        {
            add_triplets(triplets, hym_triplets, i * n_bands_sc, iy * n_bands_sc);
        }

        iy = i + 1;
        if (iy % n_ky != 0)
        {
            add_triplets(triplets, hyp_triplets, i * n_bands_sc, iy * n_bands_sc);
        }

        ix = i - n_ky;
        if (ix >= 0)
        {
            add_triplets(triplets, hxm_triplets, i * n_bands_sc, ix * n_bands_sc);
        }

        ix = i + n_ky;
        if (ix < n_kx * n_ky)
        {
            add_triplets(triplets, hxp_triplets, i * n_bands_sc, ix * n_bands_sc);
        }
    }

    return triplets;
}