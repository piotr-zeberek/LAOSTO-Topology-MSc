#ifndef SYSTEM2DELEMENTBASED_H
#define SYSTEM2DELEMENTBASED_H

#include "System2D.h"

#include <functional>

using ElementFunction = std::function<std::complex<double>(double, double)>;
using HamiltonianFunction = std::function<Hamiltonian(double, double)>;

struct System2DElementBased : public System2D
{
    using HamiltonianFunctionPointer = Hamiltonian (System2DElementBased::*)(double, double) const;
    using TripletsFunctionPointer = Triplets (System2DElementBased::*)(double, double) const;

    enum class MatrixType
    {
        Hk,
        Delta,
        Delta_Adjoint,
        mHmkT,
        HBdG
    };

    struct Element
    {
        std::size_t i;
        std::size_t j;
        ElementFunction f;

        Element() = delete;
        Element(std::size_t i, std::size_t j, ElementFunction f = [](double, double)
                                              { return 0.0; }) : i(i), j(j), f(f) {}

        std::complex<double> operator()(double x, double y) const
        {
            return f(x, y);
        }
    };

    virtual void set_Hk_elements() = 0;
    virtual void set_Delta_elements() = 0;
    virtual void set_Delta_Adjoint_elements() = 0;
    virtual void set_mHmkT_elements() = 0;

    virtual void set_Hk_discrete_ky_elements() = 0;
    virtual void set_Delta_discrete_ky_elements() = 0;
    virtual void set_Delta_Adjoint_discrete_ky_elements() = 0;
    virtual void set_mHmkT_discrete_ky_elements() = 0;

    virtual void set_Hk_discrete_elements() = 0;
    virtual void set_Delta_discrete_elements() = 0;
    virtual void set_Delta_Adjoint_discrete_elements() = 0;
    virtual void set_mHmkT_discrete_elements() = 0;

    // Matrices
    // continues hamiltonians
    Hamiltonian Hk(double kx, double ky) const;
    Hamiltonian Delta(double kx, double ky) const;
    Hamiltonian Delta_Adjoint(double kx, double ky) const;
    Hamiltonian mHmkT(double kx, double ky) const;
    Hamiltonian HBdG(double kx, double ky) const;

    // continues in kx, discretized in ky hamiltonians
    Hamiltonian Hk_discrete_ky_onsite(double kx, double y) const;
    Hamiltonian Hk_discrete_ky_hopping_p(double kx, double y) const;
    Hamiltonian Hk_discrete_ky_hopping_m(double kx, double y) const;
    Hamiltonian Hk_discrete_ky(double kx, std::size_t n_ky) const;

    Hamiltonian Delta_discrete_ky_onsite(double kx, double y) const;
    Hamiltonian Delta_discrete_ky_hopping_p(double kx, double y) const;
    Hamiltonian Delta_discrete_ky_hopping_m(double kx, double y) const;
    Hamiltonian Delta_discrete_ky(double kx, std::size_t n_ky) const;

    Hamiltonian Delta_Adjoint_discrete_ky_onsite(double kx, double y) const;
    Hamiltonian Delta_Adjoint_discrete_ky_hopping_p(double kx, double y) const;
    Hamiltonian Delta_Adjoint_discrete_ky_hopping_m(double kx, double y) const;
    Hamiltonian Delta_Adjoint_discrete_ky(double kx, std::size_t n_ky) const;

    Hamiltonian mHmkT_discrete_ky_onsite(double kx, double y) const;
    Hamiltonian mHmkT_discrete_ky_hopping_p(double kx, double y) const;
    Hamiltonian mHmkT_discrete_ky_hopping_m(double kx, double y) const;
    Hamiltonian mHmkT_discrete_ky(double kx, std::size_t n_ky) const;

    Hamiltonian HBdG_discrete_ky_onsite(double kx, double y) const;
    Hamiltonian HBdG_discrete_ky_hopping_p(double kx, double y) const;
    Hamiltonian HBdG_discrete_ky_hopping_m(double kx, double y) const;
    Hamiltonian HBdG_discrete_ky(double kx, std::size_t n_ky) const;

    // discretized in kx, discretized in ky hamiltonians
    Hamiltonian Hk_discrete_onsite(double x, double y) const;
    Hamiltonian Hk_discrete_hopping_xp(double x, double y) const;
    Hamiltonian Hk_discrete_hopping_xm(double x, double y) const;
    Hamiltonian Hk_discrete_hopping_yp(double x, double y) const;
    Hamiltonian Hk_discrete_hopping_ym(double x, double y) const;
    Hamiltonian Hk_discrete(std::size_t n_kx, std::size_t n_ky) const;

    Hamiltonian Delta_discrete_onsite(double x, double y) const;
    Hamiltonian Delta_discrete_hopping_xp(double x, double y) const;
    Hamiltonian Delta_discrete_hopping_xm(double x, double y) const;
    Hamiltonian Delta_discrete_hopping_yp(double x, double y) const;
    Hamiltonian Delta_discrete_hopping_ym(double x, double y) const;
    Hamiltonian Delta_discrete(std::size_t n_kx, std::size_t n_ky) const;

    Hamiltonian Delta_Adjoint_discrete_onsite(double x, double y) const;
    Hamiltonian Delta_Adjoint_discrete_hopping_xp(double x, double y) const;
    Hamiltonian Delta_Adjoint_discrete_hopping_xm(double x, double y) const;
    Hamiltonian Delta_Adjoint_discrete_hopping_yp(double x, double y) const;
    Hamiltonian Delta_Adjoint_discrete_hopping_ym(double x, double y) const;
    Hamiltonian Delta_Adjoint_discrete(std::size_t n_kx, std::size_t n_ky) const;

    Hamiltonian mHmkT_discrete_onsite(double x, double y) const;
    Hamiltonian mHmkT_discrete_hopping_xp(double x, double y) const;
    Hamiltonian mHmkT_discrete_hopping_xm(double x, double y) const;
    Hamiltonian mHmkT_discrete_hopping_yp(double x, double y) const;
    Hamiltonian mHmkT_discrete_hopping_ym(double x, double y) const;
    Hamiltonian mHmkT_discrete(std::size_t n_kx, std::size_t n_ky) const;

    Hamiltonian HBdG_discrete_onsite(double x, double y) const;
    Hamiltonian HBdG_discrete_hopping_xp(double x, double y) const;
    Hamiltonian HBdG_discrete_hopping_xm(double x, double y) const;
    Hamiltonian HBdG_discrete_hopping_yp(double x, double y) const;
    Hamiltonian HBdG_discrete_hopping_ym(double x, double y) const;
    Hamiltonian HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const;

    // Triplets
    // continues hamiltonians
    Triplets triplets_Hk(double kx, double ky) const;
    Triplets triplets_Delta(double kx, double ky) const;
    Triplets triplets_Delta_Adjoint(double kx, double ky) const;
    Triplets triplets_mHmkT(double kx, double ky) const;
    Triplets triplets_HBdG(double kx, double ky) const;

    // continues in kx, discretized in ky hamiltonians
    Triplets triplets_Hk_discrete_ky_onsite(double kx, double y) const;
    Triplets triplets_Hk_discrete_ky_hopping_p(double kx, double y) const;
    Triplets triplets_Hk_discrete_ky_hopping_m(double kx, double y) const;
    Triplets triplets_Hk_discrete_ky(double kx, std::size_t n_ky) const;

    Triplets triplets_Delta_discrete_ky_onsite(double kx, double y) const;
    Triplets triplets_Delta_discrete_ky_hopping_p(double kx, double y) const;
    Triplets triplets_Delta_discrete_ky_hopping_m(double kx, double y) const;
    Triplets triplets_Delta_discrete_ky(double kx, std::size_t n_ky) const;

    Triplets triplets_Delta_Adjoint_discrete_ky_onsite(double kx, double y) const;
    Triplets triplets_Delta_Adjoint_discrete_ky_hopping_p(double kx, double y) const;
    Triplets triplets_Delta_Adjoint_discrete_ky_hopping_m(double kx, double y) const;
    Triplets triplets_Delta_Adjoint_discrete_ky(double kx, std::size_t n_ky) const;

    Triplets triplets_mHmkT_discrete_ky_onsite(double kx, double y) const;
    Triplets triplets_mHmkT_discrete_ky_hopping_p(double kx, double y) const;
    Triplets triplets_mHmkT_discrete_ky_hopping_m(double kx, double y) const;
    Triplets triplets_mHmkT_discrete_ky(double kx, std::size_t n_ky) const;

    Triplets triplets_HBdG_discrete_ky_onsite(double kx, double y) const;
    Triplets triplets_HBdG_discrete_ky_hopping_p(double kx, double y) const;
    Triplets triplets_HBdG_discrete_ky_hopping_m(double kx, double y) const;
    Triplets triplets_HBdG_discrete_ky(double kx, std::size_t n_ky) const;

    // discretized in kx, discretized in ky hamiltonians
    Triplets triplets_Hk_discrete_onsite(double x, double y) const;
    Triplets triplets_Hk_discrete_hopping_xp(double x, double y) const;
    Triplets triplets_Hk_discrete_hopping_xm(double x, double y) const;
    Triplets triplets_Hk_discrete_hopping_yp(double x, double y) const;
    Triplets triplets_Hk_discrete_hopping_ym(double x, double y) const;
    Triplets triplets_Hk_discrete(std::size_t n_kx, std::size_t n_ky) const;

    Triplets triplets_Delta_discrete_onsite(double x, double y) const;
    Triplets triplets_Delta_discrete_hopping_xp(double x, double y) const;
    Triplets triplets_Delta_discrete_hopping_xm(double x, double y) const;
    Triplets triplets_Delta_discrete_hopping_yp(double x, double y) const;
    Triplets triplets_Delta_discrete_hopping_ym(double x, double y) const;
    Triplets triplets_Delta_discrete(std::size_t n_kx, std::size_t n_ky) const;

    Triplets triplets_Delta_Adjoint_discrete_onsite(double x, double y) const;
    Triplets triplets_Delta_Adjoint_discrete_hopping_xp(double x, double y) const;
    Triplets triplets_Delta_Adjoint_discrete_hopping_xm(double x, double y) const;
    Triplets triplets_Delta_Adjoint_discrete_hopping_yp(double x, double y) const;
    Triplets triplets_Delta_Adjoint_discrete_hopping_ym(double x, double y) const;
    Triplets triplets_Delta_Adjoint_discrete(std::size_t n_kx, std::size_t n_ky) const;

    Triplets triplets_mHmkT_discrete_onsite(double x, double y) const;
    Triplets triplets_mHmkT_discrete_hopping_xp(double x, double y) const;
    Triplets triplets_mHmkT_discrete_hopping_xm(double x, double y) const;
    Triplets triplets_mHmkT_discrete_hopping_yp(double x, double y) const;
    Triplets triplets_mHmkT_discrete_hopping_ym(double x, double y) const;
    Triplets triplets_mHmkT_discrete(std::size_t n_kx, std::size_t n_ky) const;

    Triplets triplets_HBdG_discrete_onsite(double x, double y) const;
    Triplets triplets_HBdG_discrete_hopping_xp(double x, double y) const;
    Triplets triplets_HBdG_discrete_hopping_xm(double x, double y) const;
    Triplets triplets_HBdG_discrete_hopping_yp(double x, double y) const;
    Triplets triplets_HBdG_discrete_hopping_ym(double x, double y) const;
    Triplets triplets_HBdG_discrete(std::size_t n_kx, std::size_t n_ky) const;

protected:
    std::vector<Element> Hk_elements;
    std::vector<Element> Delta_elements;
    std::vector<Element> Delta_Adjoint_elements;
    std::vector<Element> mHmkT_elements;

    std::vector<Element> Hk_discrete_ky_onsite_elements;
    std::vector<Element> Hk_discrete_ky_hopping_p_elements;
    std::vector<Element> Hk_discrete_ky_hopping_m_elements;

    std::vector<Element> Delta_discrete_ky_onsite_elements;
    std::vector<Element> Delta_discrete_ky_hopping_p_elements;
    std::vector<Element> Delta_discrete_ky_hopping_m_elements;

    std::vector<Element> Delta_Adjoint_discrete_ky_onsite_elements;
    std::vector<Element> Delta_Adjoint_discrete_ky_hopping_p_elements;
    std::vector<Element> Delta_Adjoint_discrete_ky_hopping_m_elements;

    std::vector<Element> mHmkT_discrete_ky_onsite_elements;
    std::vector<Element> mHmkT_discrete_ky_hopping_p_elements;
    std::vector<Element> mHmkT_discrete_ky_hopping_m_elements;

    std::vector<Element> Hk_discrete_onsite_elements;
    std::vector<Element> Hk_discrete_hopping_xp_elements;
    std::vector<Element> Hk_discrete_hopping_xm_elements;
    std::vector<Element> Hk_discrete_hopping_yp_elements;
    std::vector<Element> Hk_discrete_hopping_ym_elements;

    std::vector<Element> Delta_discrete_onsite_elements;
    std::vector<Element> Delta_discrete_hopping_xp_elements;
    std::vector<Element> Delta_discrete_hopping_xm_elements;
    std::vector<Element> Delta_discrete_hopping_yp_elements;
    std::vector<Element> Delta_discrete_hopping_ym_elements;

    std::vector<Element> Delta_Adjoint_discrete_onsite_elements;
    std::vector<Element> Delta_Adjoint_discrete_hopping_xp_elements;
    std::vector<Element> Delta_Adjoint_discrete_hopping_xm_elements;
    std::vector<Element> Delta_Adjoint_discrete_hopping_yp_elements;
    std::vector<Element> Delta_Adjoint_discrete_hopping_ym_elements;

    std::vector<Element> mHmkT_discrete_onsite_elements;
    std::vector<Element> mHmkT_discrete_hopping_xp_elements;
    std::vector<Element> mHmkT_discrete_hopping_xm_elements;
    std::vector<Element> mHmkT_discrete_hopping_yp_elements;
    std::vector<Element> mHmkT_discrete_hopping_ym_elements;

private:
    Hamiltonian assemble_matrix_from_elements(const std::vector<Element> &elements, std::size_t n_bands, double kx, double ky) const;
    Triplets assemble_triplets_from_elements(const std::vector<Element> &elements, double kx, double ky) const;

    Hamiltonian assemble_HBdG(const Hamiltonian &Hk_mat, const Hamiltonian &Delta_mat, const Hamiltonian &Delta_Adjoint_mat, const Hamiltonian &mHmkT_mat) const;
    Triplets assemble_triplets_HBdG(const Triplets &Hk_triplets, const Triplets &Delta_triplets, const Triplets &Delta_adjoint_triplets, const Triplets &mHmkT_triplets) const;

    virtual Hamiltonian assemble_matrix_discrete_ky(double kx, std::size_t n_ky, MatrixType matrix_type) const;
    virtual Hamiltonian assemble_matrix_discrete(std::size_t n_kx, std::size_t n_ky, MatrixType matrix_type) const;

    virtual Triplets assemble_triplets_discrete_ky(double kx, std::size_t n_ky, MatrixType matrix_type) const;
    virtual Triplets assemble_triplets_discrete(std::size_t n_kx, std::size_t n_ky, MatrixType matrix_type) const;
};

#endif