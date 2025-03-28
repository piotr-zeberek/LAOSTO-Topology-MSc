#include "System2DElementBased.h"
#include "utils.h"

Hamiltonian System2DElementBased::Hk(double kx, double ky) const
{
    return assemble_matrix_from_elements(Hk_elements, n_bands, kx, ky);
}

Hamiltonian System2DElementBased::Delta(double kx, double ky) const
{
    return assemble_matrix_from_elements(Delta_elements, n_bands, kx, ky);
}

Hamiltonian System2DElementBased::Delta_Adjoint(double kx, double ky) const
{
    return assemble_matrix_from_elements(Delta_Adjoint_elements, n_bands, kx, ky);
}

Hamiltonian System2DElementBased::mHmkT(double kx, double ky) const
{
    return assemble_matrix_from_elements(mHmkT_elements, n_bands, kx, ky);
}

Hamiltonian System2DElementBased::HBdG(double kx, double ky) const
{
    return assemble_HBdG(Hk(kx, ky), Delta(kx, ky), Delta_Adjoint(kx, ky), mHmkT(kx, ky));
}

Hamiltonian System2DElementBased::Hk_discrete_ky_onsite(double kx, double y) const
{
    return assemble_matrix_from_elements(Hk_discrete_ky_onsite_elements, n_bands, kx, y);
}

Hamiltonian System2DElementBased::Hk_discrete_ky_hopping_p(double kx, double y) const
{
    return assemble_matrix_from_elements(Hk_discrete_ky_hopping_p_elements, n_bands, kx, y);
}

Hamiltonian System2DElementBased::Hk_discrete_ky_hopping_m(double kx, double y) const
{
    return assemble_matrix_from_elements(Hk_discrete_ky_hopping_m_elements, n_bands, kx, y);
}

Hamiltonian System2DElementBased::Hk_discrete_ky(double kx, std::size_t n_ky) const
{
    return assemble_matrix_discrete_ky(kx, n_ky, MatrixType::Hk);
}

Hamiltonian System2DElementBased::Delta_discrete_ky_onsite(double kx, double y) const
{
    return assemble_matrix_from_elements(Delta_discrete_ky_onsite_elements, n_bands, kx, y);
}

Hamiltonian System2DElementBased::Delta_discrete_ky_hopping_p(double kx, double y) const
{
    return assemble_matrix_from_elements(Delta_discrete_ky_hopping_p_elements, n_bands, kx, y);
}

Hamiltonian System2DElementBased::Delta_discrete_ky_hopping_m(double kx, double y) const
{
    return assemble_matrix_from_elements(Delta_discrete_ky_hopping_m_elements, n_bands, kx, y);
}

Hamiltonian System2DElementBased::Delta_discrete_ky(double kx, std::size_t n_ky) const
{
    return assemble_matrix_discrete_ky(kx, n_ky, MatrixType::Delta);
}

Hamiltonian System2DElementBased::Delta_Adjoint_discrete_ky_onsite(double kx, double y) const
{
    return assemble_matrix_from_elements(Delta_Adjoint_discrete_ky_onsite_elements, n_bands, kx, y);
}

Hamiltonian System2DElementBased::Delta_Adjoint_discrete_ky_hopping_p(double kx, double y) const
{
    return assemble_matrix_from_elements(Delta_Adjoint_discrete_ky_hopping_p_elements, n_bands, kx, y);
}

Hamiltonian System2DElementBased::Delta_Adjoint_discrete_ky_hopping_m(double kx, double y) const
{
    return assemble_matrix_from_elements(Delta_Adjoint_discrete_ky_hopping_m_elements, n_bands, kx, y);
}

Hamiltonian System2DElementBased::Delta_Adjoint_discrete_ky(double kx, std::size_t n_ky) const
{
    return assemble_matrix_discrete_ky(kx, n_ky, MatrixType::Delta_Adjoint);
}

Hamiltonian System2DElementBased::mHmkT_discrete_ky_onsite(double kx, double y) const
{
    return assemble_matrix_from_elements(mHmkT_discrete_ky_onsite_elements, n_bands, kx, y);
}

Hamiltonian System2DElementBased::mHmkT_discrete_ky_hopping_p(double kx, double y) const
{
    return assemble_matrix_from_elements(mHmkT_discrete_ky_hopping_p_elements, n_bands, kx, y);
}

Hamiltonian System2DElementBased::mHmkT_discrete_ky_hopping_m(double kx, double y) const
{
    return assemble_matrix_from_elements(mHmkT_discrete_ky_hopping_m_elements, n_bands, kx, y);
}

Hamiltonian System2DElementBased::mHmkT_discrete_ky(double kx, std::size_t n_ky) const
{
    return assemble_matrix_discrete_ky(kx, n_ky, MatrixType::mHmkT);
}

Hamiltonian System2DElementBased::HBdG_discrete_ky_onsite(double kx, double y) const
{
    return assemble_HBdG(Hk_discrete_ky_onsite(kx, y), Delta_discrete_ky_onsite(kx, y), Delta_Adjoint_discrete_ky_onsite(kx, y), mHmkT_discrete_ky_onsite(kx, y));
}

Hamiltonian System2DElementBased::HBdG_discrete_ky_hopping_p(double kx, double y) const
{
    return assemble_HBdG(Hk_discrete_ky_hopping_p(kx, y), Delta_discrete_ky_hopping_p(kx, y), Delta_Adjoint_discrete_ky_hopping_p(kx, y), mHmkT_discrete_ky_hopping_p(kx, y));
}

Hamiltonian System2DElementBased::HBdG_discrete_ky_hopping_m(double kx, double y) const
{
    return assemble_HBdG(Hk_discrete_ky_hopping_m(kx, y), Delta_discrete_ky_hopping_m(kx, y), Delta_Adjoint_discrete_ky_hopping_m(kx, y), mHmkT_discrete_ky_hopping_m(kx, y));
}

Hamiltonian System2DElementBased::HBdG_discrete_ky(double kx, std::size_t n_ky) const
{
    return assemble_matrix_discrete_ky(kx, n_ky, MatrixType::HBdG);
}

Hamiltonian System2DElementBased::Hk_discrete_onsite(double x, double y) const
{
    return assemble_matrix_from_elements(Hk_discrete_onsite_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::Hk_discrete_hopping_xp(double x, double y) const
{
    return assemble_matrix_from_elements(Hk_discrete_hopping_xp_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::Hk_discrete_hopping_xm(double x, double y) const
{
    return assemble_matrix_from_elements(Hk_discrete_hopping_xm_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::Hk_discrete_hopping_yp(double x, double y) const
{
    return assemble_matrix_from_elements(Hk_discrete_hopping_yp_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::Hk_discrete_hopping_ym(double x, double y) const
{
    return assemble_matrix_from_elements(Hk_discrete_hopping_ym_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::Hk_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_matrix_discrete(n_kx, n_ky, MatrixType::Hk);
}

Hamiltonian System2DElementBased::Delta_discrete_onsite(double x, double y) const
{
    return assemble_matrix_from_elements(Delta_discrete_onsite_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::Delta_discrete_hopping_xp(double x, double y) const
{
    return assemble_matrix_from_elements(Delta_discrete_hopping_xp_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::Delta_discrete_hopping_xm(double x, double y) const
{
    return assemble_matrix_from_elements(Delta_discrete_hopping_xm_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::Delta_discrete_hopping_yp(double x, double y) const
{
    return assemble_matrix_from_elements(Delta_discrete_hopping_yp_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::Delta_discrete_hopping_ym(double x, double y) const
{
    return assemble_matrix_from_elements(Delta_discrete_hopping_ym_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::Delta_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_matrix_discrete(n_kx, n_ky, MatrixType::Delta);
}

Hamiltonian System2DElementBased::Delta_Adjoint_discrete_onsite(double x, double y) const
{
    return assemble_matrix_from_elements(Delta_Adjoint_discrete_onsite_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::Delta_Adjoint_discrete_hopping_xp(double x, double y) const
{
    return assemble_matrix_from_elements(Delta_Adjoint_discrete_hopping_xp_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::Delta_Adjoint_discrete_hopping_xm(double x, double y) const
{
    return assemble_matrix_from_elements(Delta_Adjoint_discrete_hopping_xm_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::Delta_Adjoint_discrete_hopping_yp(double x, double y) const
{
    return assemble_matrix_from_elements(Delta_Adjoint_discrete_hopping_yp_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::Delta_Adjoint_discrete_hopping_ym(double x, double y) const
{
    return assemble_matrix_from_elements(Delta_Adjoint_discrete_hopping_ym_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::Delta_Adjoint_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_matrix_discrete(n_kx, n_ky, MatrixType::Delta_Adjoint);
}

Hamiltonian System2DElementBased::mHmkT_discrete_onsite(double x, double y) const
{
    return assemble_matrix_from_elements(mHmkT_discrete_onsite_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::mHmkT_discrete_hopping_xp(double x, double y) const
{
    return assemble_matrix_from_elements(mHmkT_discrete_hopping_xp_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::mHmkT_discrete_hopping_xm(double x, double y) const
{
    return assemble_matrix_from_elements(mHmkT_discrete_hopping_xm_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::mHmkT_discrete_hopping_yp(double x, double y) const
{
    return assemble_matrix_from_elements(mHmkT_discrete_hopping_yp_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::mHmkT_discrete_hopping_ym(double x, double y) const
{
    return assemble_matrix_from_elements(mHmkT_discrete_hopping_ym_elements, n_bands, x, y);
}

Hamiltonian System2DElementBased::mHmkT_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_matrix_discrete(n_kx, n_ky, MatrixType::mHmkT);
}

Hamiltonian System2DElementBased::HBdG_discrete_onsite(double x, double y) const
{
    return assemble_HBdG(Hk_discrete_onsite(x, y), Delta_discrete_onsite(x, y), Delta_Adjoint_discrete_onsite(x, y), mHmkT_discrete_onsite(x, y));
}

Hamiltonian System2DElementBased::HBdG_discrete_hopping_xp(double x, double y) const
{
    return assemble_HBdG(Hk_discrete_hopping_xp(x, y), Delta_discrete_hopping_xp(x, y), Delta_Adjoint_discrete_hopping_xp(x, y), mHmkT_discrete_hopping_xp(x, y));
}

Hamiltonian System2DElementBased::HBdG_discrete_hopping_xm(double x, double y) const
{
    return assemble_HBdG(Hk_discrete_hopping_xm(x, y), Delta_discrete_hopping_xm(x, y), Delta_Adjoint_discrete_hopping_xm(x, y), mHmkT_discrete_hopping_xm(x, y));
}

Hamiltonian System2DElementBased::HBdG_discrete_hopping_yp(double x, double y) const
{
    return assemble_HBdG(Hk_discrete_hopping_yp(x, y), Delta_discrete_hopping_yp(x, y), Delta_Adjoint_discrete_hopping_yp(x, y), mHmkT_discrete_hopping_yp(x, y));
}

Hamiltonian System2DElementBased::HBdG_discrete_hopping_ym(double x, double y) const
{
    return assemble_HBdG(Hk_discrete_hopping_ym(x, y), Delta_discrete_hopping_ym(x, y), Delta_Adjoint_discrete_hopping_ym(x, y), mHmkT_discrete_hopping_ym(x, y));
}

Hamiltonian System2DElementBased::HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_matrix_discrete(n_kx, n_ky, MatrixType::HBdG);
}

Triplets System2DElementBased::triplets_Hk(double kx, double ky) const
{
    return assemble_triplets_from_elements(Hk_elements, kx, ky);
}

Triplets System2DElementBased::triplets_Delta(double kx, double ky) const
{
    return assemble_triplets_from_elements(Delta_elements, kx, ky);
}

Triplets System2DElementBased::triplets_Delta_Adjoint(double kx, double ky) const
{
    return assemble_triplets_from_elements(Delta_Adjoint_elements, kx, ky);
}

Triplets System2DElementBased::triplets_mHmkT(double kx, double ky) const
{
    return assemble_triplets_from_elements(mHmkT_elements, kx, ky);
}

Triplets System2DElementBased::triplets_HBdG(double kx, double ky) const
{
    return assemble_triplets_HBdG(triplets_Hk(kx, ky), triplets_Delta(kx, ky), triplets_Delta_Adjoint(kx, ky), triplets_mHmkT(kx, ky));
}

Triplets System2DElementBased::triplets_Hk_discrete_ky_onsite(double kx, double y) const
{
    return assemble_triplets_from_elements(Hk_discrete_ky_onsite_elements, kx, y);
}

Triplets System2DElementBased::triplets_Hk_discrete_ky_hopping_p(double kx, double y) const
{
    return assemble_triplets_from_elements(Hk_discrete_ky_hopping_p_elements, kx, y);
}

Triplets System2DElementBased::triplets_Hk_discrete_ky_hopping_m(double kx, double y) const
{
    return assemble_triplets_from_elements(Hk_discrete_ky_hopping_m_elements, kx, y);
}

Triplets System2DElementBased::triplets_Hk_discrete_ky(double kx, std::size_t n_ky) const
{
    return assemble_triplets_discrete_ky(kx, n_ky, MatrixType::Hk);
}

Triplets System2DElementBased::triplets_Delta_discrete_ky_onsite(double kx, double y) const
{
    return assemble_triplets_from_elements(Delta_discrete_ky_onsite_elements, kx, y);
}

Triplets System2DElementBased::triplets_Delta_discrete_ky_hopping_p(double kx, double y) const
{
    return assemble_triplets_from_elements(Delta_discrete_ky_hopping_p_elements, kx, y);
}

Triplets System2DElementBased::triplets_Delta_discrete_ky_hopping_m(double kx, double y) const
{
    return assemble_triplets_from_elements(Delta_discrete_ky_hopping_m_elements, kx, y);
}

Triplets System2DElementBased::triplets_Delta_discrete_ky(double kx, std::size_t n_ky) const
{
    return assemble_triplets_discrete_ky(kx, n_ky, MatrixType::Delta);
}

Triplets System2DElementBased::triplets_Delta_Adjoint_discrete_ky_onsite(double kx, double y) const
{
    return assemble_triplets_from_elements(Delta_Adjoint_discrete_ky_onsite_elements, kx, y);
}

Triplets System2DElementBased::triplets_Delta_Adjoint_discrete_ky_hopping_p(double kx, double y) const
{
    return assemble_triplets_from_elements(Delta_Adjoint_discrete_ky_hopping_p_elements, kx, y);
}

Triplets System2DElementBased::triplets_Delta_Adjoint_discrete_ky_hopping_m(double kx, double y) const
{
    return assemble_triplets_from_elements(Delta_Adjoint_discrete_ky_hopping_m_elements, kx, y);
}

Triplets System2DElementBased::triplets_Delta_Adjoint_discrete_ky(double kx, std::size_t n_ky) const
{
    return assemble_triplets_discrete_ky(kx, n_ky, MatrixType::Delta_Adjoint);
}

Triplets System2DElementBased::triplets_mHmkT_discrete_ky_onsite(double kx, double y) const
{
    return assemble_triplets_from_elements(mHmkT_discrete_ky_onsite_elements, kx, y);
}

Triplets System2DElementBased::triplets_mHmkT_discrete_ky_hopping_p(double kx, double y) const
{
    return assemble_triplets_from_elements(mHmkT_discrete_ky_hopping_p_elements, kx, y);
}

Triplets System2DElementBased::triplets_mHmkT_discrete_ky_hopping_m(double kx, double y) const
{
    return assemble_triplets_from_elements(mHmkT_discrete_ky_hopping_m_elements, kx, y);
}

Triplets System2DElementBased::triplets_mHmkT_discrete_ky(double kx, std::size_t n_ky) const
{
    return assemble_triplets_discrete_ky(kx, n_ky, MatrixType::mHmkT);
}

Triplets System2DElementBased::triplets_HBdG_discrete_ky_onsite(double kx, double y) const
{
    return assemble_triplets_HBdG(triplets_Hk_discrete_ky_onsite(kx, y), triplets_Delta_discrete_ky_onsite(kx, y), triplets_Delta_Adjoint_discrete_ky_onsite(kx, y), triplets_mHmkT_discrete_ky_onsite(kx, y));
}

Triplets System2DElementBased::triplets_HBdG_discrete_ky_hopping_p(double kx, double y) const
{
    return assemble_triplets_HBdG(triplets_Hk_discrete_ky_hopping_p(kx, y), triplets_Delta_discrete_ky_hopping_p(kx, y), triplets_Delta_Adjoint_discrete_ky_hopping_p(kx, y), triplets_mHmkT_discrete_ky_hopping_p(kx, y));
}

Triplets System2DElementBased::triplets_HBdG_discrete_ky_hopping_m(double kx, double y) const
{
    return assemble_triplets_HBdG(triplets_Hk_discrete_ky_hopping_m(kx, y), triplets_Delta_discrete_ky_hopping_m(kx, y), triplets_Delta_Adjoint_discrete_ky_hopping_m(kx, y), triplets_mHmkT_discrete_ky_hopping_m(kx, y));
}

Triplets System2DElementBased::triplets_HBdG_discrete_ky(double kx, std::size_t n_ky) const
{
    return assemble_triplets_discrete_ky(kx, n_ky, MatrixType::HBdG);
}

Triplets System2DElementBased::triplets_Hk_discrete_onsite(double x, double y) const
{
    return assemble_triplets_from_elements(Hk_discrete_onsite_elements, x, y);
}

Triplets System2DElementBased::triplets_Hk_discrete_hopping_xp(double x, double y) const
{
    return assemble_triplets_from_elements(Hk_discrete_hopping_xp_elements, x, y);
}

Triplets System2DElementBased::triplets_Hk_discrete_hopping_xm(double x, double y) const
{
    return assemble_triplets_from_elements(Hk_discrete_hopping_xm_elements, x, y);
}

Triplets System2DElementBased::triplets_Hk_discrete_hopping_yp(double x, double y) const
{
    return assemble_triplets_from_elements(Hk_discrete_hopping_yp_elements, x, y);
}

Triplets System2DElementBased::triplets_Hk_discrete_hopping_ym(double x, double y) const
{
    return assemble_triplets_from_elements(Hk_discrete_hopping_ym_elements, x, y);
}

Triplets System2DElementBased::triplets_Hk_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_triplets_discrete(n_kx, n_ky, MatrixType::Hk);
}

Triplets System2DElementBased::triplets_Delta_discrete_onsite(double x, double y) const
{
    return assemble_triplets_from_elements(Delta_discrete_onsite_elements, x, y);
}

Triplets System2DElementBased::triplets_Delta_discrete_hopping_xp(double x, double y) const
{
    return assemble_triplets_from_elements(Delta_discrete_hopping_xp_elements, x, y);
}

Triplets System2DElementBased::triplets_Delta_discrete_hopping_xm(double x, double y) const
{
    return assemble_triplets_from_elements(Delta_discrete_hopping_xm_elements, x, y);
}

Triplets System2DElementBased::triplets_Delta_discrete_hopping_yp(double x, double y) const
{
    return assemble_triplets_from_elements(Delta_discrete_hopping_yp_elements, x, y);
}

Triplets System2DElementBased::triplets_Delta_discrete_hopping_ym(double x, double y) const
{
    return assemble_triplets_from_elements(Delta_discrete_hopping_ym_elements, x, y);
}

Triplets System2DElementBased::triplets_Delta_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_triplets_discrete(n_kx, n_ky, MatrixType::Delta);
}

Triplets System2DElementBased::triplets_Delta_Adjoint_discrete_onsite(double x, double y) const
{
    return assemble_triplets_from_elements(Delta_Adjoint_discrete_onsite_elements, x, y);
}

Triplets System2DElementBased::triplets_Delta_Adjoint_discrete_hopping_xp(double x, double y) const
{
    return assemble_triplets_from_elements(Delta_Adjoint_discrete_hopping_xp_elements, x, y);
}

Triplets System2DElementBased::triplets_Delta_Adjoint_discrete_hopping_xm(double x, double y) const
{
    return assemble_triplets_from_elements(Delta_Adjoint_discrete_hopping_xm_elements, x, y);
}

Triplets System2DElementBased::triplets_Delta_Adjoint_discrete_hopping_yp(double x, double y) const
{
    return assemble_triplets_from_elements(Delta_Adjoint_discrete_hopping_yp_elements, x, y);
}

Triplets System2DElementBased::triplets_Delta_Adjoint_discrete_hopping_ym(double x, double y) const
{
    return assemble_triplets_from_elements(Delta_Adjoint_discrete_hopping_ym_elements, x, y);
}

Triplets System2DElementBased::triplets_Delta_Adjoint_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_triplets_discrete(n_kx, n_ky, MatrixType::Delta_Adjoint);
}

Triplets System2DElementBased::triplets_mHmkT_discrete_onsite(double x, double y) const
{
    return assemble_triplets_from_elements(mHmkT_discrete_onsite_elements, x, y);
}

Triplets System2DElementBased::triplets_mHmkT_discrete_hopping_xp(double x, double y) const
{
    return assemble_triplets_from_elements(mHmkT_discrete_hopping_xp_elements, x, y);
}

Triplets System2DElementBased::triplets_mHmkT_discrete_hopping_xm(double x, double y) const
{
    return assemble_triplets_from_elements(mHmkT_discrete_hopping_xm_elements, x, y);
}

Triplets System2DElementBased::triplets_mHmkT_discrete_hopping_yp(double x, double y) const
{
    return assemble_triplets_from_elements(mHmkT_discrete_hopping_yp_elements, x, y);
}

Triplets System2DElementBased::triplets_mHmkT_discrete_hopping_ym(double x, double y) const
{
    return assemble_triplets_from_elements(mHmkT_discrete_hopping_ym_elements, x, y);
}

Triplets System2DElementBased::triplets_mHmkT_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_triplets_discrete(n_kx, n_ky, MatrixType::mHmkT);
}

Triplets System2DElementBased::triplets_HBdG_discrete_onsite(double x, double y) const
{
    return assemble_triplets_HBdG(triplets_Hk_discrete_onsite(x, y), triplets_Delta_discrete_onsite(x, y), triplets_Delta_Adjoint_discrete_onsite(x, y), triplets_mHmkT_discrete_onsite(x, y));
}

Triplets System2DElementBased::triplets_HBdG_discrete_hopping_xp(double x, double y) const
{
    return assemble_triplets_HBdG(triplets_Hk_discrete_hopping_xp(x, y), triplets_Delta_discrete_hopping_xp(x, y), triplets_Delta_Adjoint_discrete_hopping_xp(x, y), triplets_mHmkT_discrete_hopping_xp(x, y));
}

Triplets System2DElementBased::triplets_HBdG_discrete_hopping_xm(double x, double y) const
{
    return assemble_triplets_HBdG(triplets_Hk_discrete_hopping_xm(x, y), triplets_Delta_discrete_hopping_xm(x, y), triplets_Delta_Adjoint_discrete_hopping_xm(x, y), triplets_mHmkT_discrete_hopping_xm(x, y));
}

Triplets System2DElementBased::triplets_HBdG_discrete_hopping_yp(double x, double y) const
{
    return assemble_triplets_HBdG(triplets_Hk_discrete_hopping_yp(x, y), triplets_Delta_discrete_hopping_yp(x, y), triplets_Delta_Adjoint_discrete_hopping_yp(x, y), triplets_mHmkT_discrete_hopping_yp(x, y));
}

Triplets System2DElementBased::triplets_HBdG_discrete_hopping_ym(double x, double y) const
{
    return assemble_triplets_HBdG(triplets_Hk_discrete_hopping_ym(x, y), triplets_Delta_discrete_hopping_ym(x, y), triplets_Delta_Adjoint_discrete_hopping_ym(x, y), triplets_mHmkT_discrete_hopping_ym(x, y));
}

Triplets System2DElementBased::triplets_HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const
{
    return assemble_triplets_discrete(n_kx, n_ky, MatrixType::HBdG);
}

Hamiltonian System2DElementBased::assemble_matrix_from_elements(const std::vector<Element> &elements, std::size_t n_bands, double kx, double ky) const
{
    Hamiltonian H = Hamiltonian::Zero(n_bands, n_bands);

    for (const auto &element : elements)
    {
        H(element.i, element.j) = element(kx, ky);
    }

    return H;
}

Triplets System2DElementBased::assemble_triplets_from_elements(const std::vector<Element> &elements, double kx, double ky) const
{
    Triplets triplets;
    triplets.reserve(elements.size());

    for (const auto &element : elements)
    {
        triplets.emplace_back(element.i, element.j, element(kx, ky));
    }

    return triplets;
}

Hamiltonian System2DElementBased::assemble_HBdG(const Hamiltonian &Hk_mat, const Hamiltonian &Delta_mat, const Hamiltonian &Delta_Adjoint_mat, const Hamiltonian &mHmkT_mat) const
{
    Hamiltonian H(2 * n_bands, 2 * n_bands);

    H << Hk_mat, Delta_mat,
        Delta_Adjoint_mat, mHmkT_mat;

    return H;
}

Triplets System2DElementBased::assemble_triplets_HBdG(const Triplets &Hk_triplets, const Triplets &Delta_triplets, const Triplets &Delta_adjoint_triplets, const Triplets &mHmkT_triplets) const
{
    Triplets triplets;
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

Hamiltonian System2DElementBased::assemble_matrix_discrete_ky(double kx, std::size_t n_ky, MatrixType matrix_type) const
{
    // choose the right onsite and hopping functions, use them to assemble the matrix, use pointers
    // to the functions to avoid the switch statement
    HamiltonianFunctionPointer onsite;
    HamiltonianFunctionPointer hopping_p;
    HamiltonianFunctionPointer hopping_m;

    switch (matrix_type)
    {
    case MatrixType::Hk:
        onsite = &System2DElementBased::Hk_discrete_ky_onsite;
        hopping_p = &System2DElementBased::Hk_discrete_ky_hopping_p;
        hopping_m = &System2DElementBased::Hk_discrete_ky_hopping_m;
        break;
    case MatrixType::Delta:
        onsite = &System2DElementBased::Delta_discrete_ky_onsite;
        hopping_p = &System2DElementBased::Delta_discrete_ky_hopping_p;
        hopping_m = &System2DElementBased::Delta_discrete_ky_hopping_m;
        break;
    case MatrixType::Delta_Adjoint:
        onsite = &System2DElementBased::Delta_Adjoint_discrete_ky_onsite;
        hopping_p = &System2DElementBased::Delta_Adjoint_discrete_ky_hopping_p;
        hopping_m = &System2DElementBased::Delta_Adjoint_discrete_ky_hopping_m;
        break;
    case MatrixType::mHmkT:
        onsite = &System2DElementBased::mHmkT_discrete_ky_onsite;
        hopping_p = &System2DElementBased::mHmkT_discrete_ky_hopping_p;
        hopping_m = &System2DElementBased::mHmkT_discrete_ky_hopping_m;
        break;
    case MatrixType::HBdG:
        onsite = &System2DElementBased::HBdG_discrete_ky_onsite;
        hopping_p = &System2DElementBased::HBdG_discrete_ky_hopping_p;
        hopping_m = &System2DElementBased::HBdG_discrete_ky_hopping_m;
        break;
    }

    int subsize = matrix_type == MatrixType::HBdG ? 2 * n_bands : n_bands;

    Hamiltonian H = Hamiltonian::Zero(subsize * n_ky, subsize * n_ky);

    int iy = 0;

    for (int i = 0; i < n_ky; ++i)
    {
        double y = dy * i;

        Hamiltonian o = std::invoke(onsite, this, kx, y);
        Hamiltonian hp = std::invoke(hopping_p, this, kx, y);
        Hamiltonian hm = std::invoke(hopping_m, this, kx, y);

        H.block(i * subsize, i * subsize, subsize, subsize) = o;

        iy = i - 1;
        if (i % n_ky != 0)
            H.block(i * subsize, iy * subsize, subsize, subsize) = hm;

        iy = i + 1;
        if (iy % n_ky != 0)
            H.block(i * subsize, iy * subsize, subsize, subsize) = hp;
    }

    return H;
}

Hamiltonian System2DElementBased::assemble_matrix_discrete(std::size_t n_kx, std::size_t n_ky, MatrixType matrix_type) const
{
    HamiltonianFunctionPointer onsite;
    HamiltonianFunctionPointer hopping_xp;
    HamiltonianFunctionPointer hopping_xm;
    HamiltonianFunctionPointer hopping_yp;
    HamiltonianFunctionPointer hopping_ym;

    switch (matrix_type)
    {
    case MatrixType::Hk:
        onsite = &System2DElementBased::Hk_discrete_onsite;
        hopping_xp = &System2DElementBased::Hk_discrete_hopping_xp;
        hopping_xm = &System2DElementBased::Hk_discrete_hopping_xm;
        hopping_yp = &System2DElementBased::Hk_discrete_hopping_yp;
        hopping_ym = &System2DElementBased::Hk_discrete_hopping_ym;
        break;
    case MatrixType::Delta:
        onsite = &System2DElementBased::Delta_discrete_onsite;
        hopping_xp = &System2DElementBased::Delta_discrete_hopping_xp;
        hopping_xm = &System2DElementBased::Delta_discrete_hopping_xm;
        hopping_yp = &System2DElementBased::Delta_discrete_hopping_yp;
        hopping_ym = &System2DElementBased::Delta_discrete_hopping_ym;
        break;
    case MatrixType::Delta_Adjoint:
        onsite = &System2DElementBased::Delta_Adjoint_discrete_onsite;
        hopping_xp = &System2DElementBased::Delta_Adjoint_discrete_hopping_xp;
        hopping_xm = &System2DElementBased::Delta_Adjoint_discrete_hopping_xm;
        hopping_yp = &System2DElementBased::Delta_Adjoint_discrete_hopping_yp;
        hopping_ym = &System2DElementBased::Delta_Adjoint_discrete_hopping_ym;
        break;
    case MatrixType::mHmkT:
        onsite = &System2DElementBased::mHmkT_discrete_onsite;
        hopping_xp = &System2DElementBased::mHmkT_discrete_hopping_xp;
        hopping_xm = &System2DElementBased::mHmkT_discrete_hopping_xm;
        hopping_yp = &System2DElementBased::mHmkT_discrete_hopping_yp;
        hopping_ym = &System2DElementBased::mHmkT_discrete_hopping_ym;
        break;
    case MatrixType::HBdG:
        onsite = &System2DElementBased::HBdG_discrete_onsite;
        hopping_xp = &System2DElementBased::HBdG_discrete_hopping_xp;
        hopping_xm = &System2DElementBased::HBdG_discrete_hopping_xm;
        hopping_yp = &System2DElementBased::HBdG_discrete_hopping_yp;
        hopping_ym = &System2DElementBased::HBdG_discrete_hopping_ym;
        break;
    }

    int subsize = matrix_type == MatrixType::HBdG ? 2 * n_bands : n_bands;

    Hamiltonian H = Hamiltonian::Zero(subsize * n_kx * n_ky, subsize * n_kx * n_ky);

    int ix = 0;
    int iy = 0;

    for (int i = 0; i < n_kx * n_ky; ++i)
    {
        double x = dx * (i / n_ky);
        double y = dy * (i % n_ky);

        Hamiltonian o = std::invoke(onsite, this, x, y);
        Hamiltonian hxp = std::invoke(hopping_xp, this, x, y);
        Hamiltonian hxm = std::invoke(hopping_xm, this, x, y);
        Hamiltonian hyp = std::invoke(hopping_yp, this, x, y);
        Hamiltonian hym = std::invoke(hopping_ym, this, x, y);

        H.block(i * subsize, i * subsize, subsize, subsize) = o;

        iy = i - 1;
        if (i % n_ky != 0)
            H.block(i * subsize, iy * subsize, subsize, subsize) = hym;

        iy = i + 1;
        if (iy % n_ky != 0)
            H.block(i * subsize, iy * subsize, subsize, subsize) = hyp;

        ix = i - n_ky;
        if (ix >= 0)
            H.block(i * subsize, ix * subsize, subsize, subsize) = hxm;

        ix = i + n_ky;
        if (ix < n_kx * n_ky)
            H.block(i * subsize, ix * subsize, subsize, subsize) = hxp;
    }

    return H;
}

Triplets System2DElementBased::assemble_triplets_discrete_ky(double kx, std::size_t n_ky, MatrixType matrix_type) const
{
    TripletsFunctionPointer onsite;
    TripletsFunctionPointer hopping_p;
    TripletsFunctionPointer hopping_m;

    switch (matrix_type)
    {
    case MatrixType::Hk:
        onsite = &System2DElementBased::triplets_Hk_discrete_ky_onsite;
        hopping_p = &System2DElementBased::triplets_Hk_discrete_ky_hopping_p;
        hopping_m = &System2DElementBased::triplets_Hk_discrete_ky_hopping_m;
        break;
    case MatrixType::Delta:
        onsite = &System2DElementBased::triplets_Delta_discrete_ky_onsite;
        hopping_p = &System2DElementBased::triplets_Delta_discrete_ky_hopping_p;
        hopping_m = &System2DElementBased::triplets_Delta_discrete_ky_hopping_m;
        break;
    case MatrixType::Delta_Adjoint:
        onsite = &System2DElementBased::triplets_Delta_Adjoint_discrete_ky_onsite;
        hopping_p = &System2DElementBased::triplets_Delta_Adjoint_discrete_ky_hopping_p;
        hopping_m = &System2DElementBased::triplets_Delta_Adjoint_discrete_ky_hopping_m;
        break;
    case MatrixType::mHmkT:
        onsite = &System2DElementBased::triplets_mHmkT_discrete_ky_onsite;
        hopping_p = &System2DElementBased::triplets_mHmkT_discrete_ky_hopping_p;
        hopping_m = &System2DElementBased::triplets_mHmkT_discrete_ky_hopping_m;
        break;
    case MatrixType::HBdG:
        onsite = &System2DElementBased::triplets_HBdG_discrete_ky_onsite;
        hopping_p = &System2DElementBased::triplets_HBdG_discrete_ky_hopping_p;
        hopping_m = &System2DElementBased::triplets_HBdG_discrete_ky_hopping_m;
        break;
    }

    int subsize = matrix_type == MatrixType::HBdG ? 2 * n_bands : n_bands;

    Triplets triplets;
    triplets.reserve(subsize * subsize * n_ky * 3);

    int iy = 0;

    for (int i = 0; i < n_ky; ++i)
    {
        double y = dy * i;

        Triplets triplets_o = std::invoke(onsite, this, kx, y);
        Triplets triplets_hp = std::invoke(hopping_p, this, kx, y);
        Triplets triplets_hm = std::invoke(hopping_m, this, kx, y);

        add_triplets(triplets, triplets_o, i * subsize, i * subsize);

        iy = i - 1;
        if (i % n_ky != 0)
            add_triplets(triplets, triplets_hm, i * subsize, iy * subsize);

        iy = i + 1;
        if (iy % n_ky != 0)
            add_triplets(triplets, triplets_hp, i * subsize, iy * subsize);
    }

    return triplets;
}

Triplets System2DElementBased::assemble_triplets_discrete(std::size_t n_kx, std::size_t n_ky, MatrixType matrix_type) const
{
    TripletsFunctionPointer onsite;
    TripletsFunctionPointer hopping_xp;
    TripletsFunctionPointer hopping_xm;
    TripletsFunctionPointer hopping_yp;
    TripletsFunctionPointer hopping_ym;

    switch (matrix_type)
    {
    case MatrixType::Hk:
        onsite = &System2DElementBased::triplets_Hk_discrete_onsite;
        hopping_xp = &System2DElementBased::triplets_Hk_discrete_hopping_xp;
        hopping_xm = &System2DElementBased::triplets_Hk_discrete_hopping_xm;
        hopping_yp = &System2DElementBased::triplets_Hk_discrete_hopping_yp;
        hopping_ym = &System2DElementBased::triplets_Hk_discrete_hopping_ym;
        break;
    case MatrixType::Delta:
        onsite = &System2DElementBased::triplets_Delta_discrete_onsite;
        hopping_xp = &System2DElementBased::triplets_Delta_discrete_hopping_xp;
        hopping_xm = &System2DElementBased::triplets_Delta_discrete_hopping_xm;
        hopping_yp = &System2DElementBased::triplets_Delta_discrete_hopping_yp;
        hopping_ym = &System2DElementBased::triplets_Delta_discrete_hopping_ym;
        break;
    case MatrixType::Delta_Adjoint:
        onsite = &System2DElementBased::triplets_Delta_Adjoint_discrete_onsite;
        hopping_xp = &System2DElementBased::triplets_Delta_Adjoint_discrete_hopping_xp;
        hopping_xm = &System2DElementBased::triplets_Delta_Adjoint_discrete_hopping_xm;
        hopping_yp = &System2DElementBased::triplets_Delta_Adjoint_discrete_hopping_yp;
        hopping_ym = &System2DElementBased::triplets_Delta_Adjoint_discrete_hopping_ym;
        break;
    case MatrixType::mHmkT:
        onsite = &System2DElementBased::triplets_mHmkT_discrete_onsite;
        hopping_xp = &System2DElementBased::triplets_mHmkT_discrete_hopping_xp;
        hopping_xm = &System2DElementBased::triplets_mHmkT_discrete_hopping_xm;
        hopping_yp = &System2DElementBased::triplets_mHmkT_discrete_hopping_yp;
        hopping_ym = &System2DElementBased::triplets_mHmkT_discrete_hopping_ym;
        break;
    case MatrixType::HBdG:
        onsite = &System2DElementBased::triplets_HBdG_discrete_onsite;
        hopping_xp = &System2DElementBased::triplets_HBdG_discrete_hopping_xp;
        hopping_xm = &System2DElementBased::triplets_HBdG_discrete_hopping_xm;
        hopping_yp = &System2DElementBased::triplets_HBdG_discrete_hopping_yp;
        hopping_ym = &System2DElementBased::triplets_HBdG_discrete_hopping_ym;
        break;
    }

    int subsize = matrix_type == MatrixType::HBdG ? 2 * n_bands : n_bands;

    Triplets triplets;
    triplets.reserve(subsize * subsize * n_kx * n_ky * 5);

    int ix = 0;
    int iy = 0;

    for (int i = 0; i < n_kx * n_ky; ++i)
    {
        double x = dx * (i / n_ky);
        double y = dy * (i % n_ky);

        Triplets triplets_o = std::invoke(onsite, this, x, y);
        Triplets triplets_hxp = std::invoke(hopping_xp, this, x, y);
        Triplets triplets_hxm = std::invoke(hopping_xm, this, x, y);
        Triplets triplets_hyp = std::invoke(hopping_yp, this, x, y);
        Triplets triplets_hym = std::invoke(hopping_ym, this, x, y);

        add_triplets(triplets, triplets_o, i * subsize, i * subsize);

        iy = i - 1;
        if (i % n_ky != 0)
            add_triplets(triplets, triplets_hym, i * subsize, iy * subsize);

        iy = i + 1;
        if (iy % n_ky != 0)
            add_triplets(triplets, triplets_hyp, i * subsize, iy * subsize);

        ix = i - n_ky;
        if (ix >= 0)
            add_triplets(triplets, triplets_hxm, i * subsize, ix * subsize);

        ix = i + n_ky;
        if (ix < n_kx * n_ky)
            add_triplets(triplets, triplets_hxp, i * subsize, ix * subsize);
    }

    return triplets;
}
