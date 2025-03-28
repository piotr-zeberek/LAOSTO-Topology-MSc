#include "ToyModelElementBased.h"

void ToyModelElementBased::set_Hk_elements()
{
    Element Hk00(0, 0, [this](double kx, double ky)
                 { return Ek(kx, ky) + HBz() - mu; });

    Element Hk01(0, 1, [this](double kx, double ky)
                 { return HBx_iy() + HRSO(kx, ky); });

    Element Hk10(1, 0, [this](double kx, double ky)
                 { return std::conj(HBx_iy()) + std::conj(HRSO(kx, ky)); });

    Element Hk11(1, 1, [this](double kx, double ky)
                 { return Ek(kx, ky) - HBz() - mu; });

    Hk_elements = {Hk00, Hk01, Hk10, Hk11};
}

void ToyModelElementBased::set_Delta_elements()
{
    Element Delta01(0, 1, [this](double kx, double ky)
                    { return delta_SC; });

    Element Delta10(1, 0, [this](double kx, double ky)
                    { return -delta_SC; });

    Delta_elements = {Delta01, Delta10};
}

void ToyModelElementBased::set_Delta_Adjoint_elements()
{
    Element Delta01(0, 1, [this](double kx, double ky)
                    { return -delta_SC; });

    Element Delta10(1, 0, [this](double kx, double ky)
                    { return delta_SC; });

    Delta_Adjoint_elements = {Delta01, Delta10};
}

void ToyModelElementBased::set_mHmkT_elements()
{
    Element Hk00(0, 0, [this](double kx, double ky)
                 { return -Ek(kx, ky) - HBz() + mu; });

    Element Hk01(0, 1, [this](double kx, double ky)
                 { return -std::conj(HBx_iy()) + std::conj(HRSO(kx, ky)); });

    Element Hk10(1, 0, [this](double kx, double ky)
                 { return -HBx_iy() + HRSO(kx, ky); });

    Element Hk11(1, 1, [this](double kx, double ky)
                 { return -Ek(kx, ky) + HBz() + mu; });

    mHmkT_elements = {Hk00, Hk01, Hk10, Hk11};
}

void ToyModelElementBased::set_Hk_discrete_ky_elements()
{

    Element onsite00(0, 0, [this](double kx, double y)
                     { return Ek(kx) + HBz() - mu; });

    Element onsite01(0, 1, [this](double kx, double y)
                     { return HBx_iy() + HRSO(kx); });

    Element onsite10(1, 0, [this](double kx, double y)
                     { return std::conj(HBx_iy()) + std::conj(HRSO(kx)); });

    Element onsite11(1, 1, [this](double kx, double y)
                     { return Ek(kx) - HBz() - mu; });

    Hk_discrete_ky_onsite_elements = {onsite00, onsite01, onsite10, onsite11};

    Element hopping_p00(0, 0, [this](double kx, double y)
                        { return -ty; });

    Element hopping_p01(0, 1, [this](double kx, double y)
                        { return -1i * 0.5 * delta_RSO_y; });

    Element hopping_p10(1, 0, [this](double kx, double y)
                        { return -1i * 0.5 * delta_RSO_y; });

    Element hopping_p11(1, 1, [this](double kx, double y)
                        { return -ty; });

    Hk_discrete_ky_hopping_p_elements = {hopping_p00, hopping_p01, hopping_p10, hopping_p11};

    Element hopping_m00(0, 0, [this](double kx, double y)
                        { return -ty; });

    Element hopping_m01(0, 1, [this](double kx, double y)
                        { return 1i * 0.5 * delta_RSO_y; });

    Element hopping_m10(1, 0, [this](double kx, double y)
                        { return 1i * 0.5 * delta_RSO_y; });

    Element hopping_m11(1, 1, [this](double kx, double y)
                        { return -ty; });

    Hk_discrete_ky_hopping_m_elements = {hopping_m00, hopping_m01, hopping_m10, hopping_m11};
}

void ToyModelElementBased::set_Delta_discrete_ky_elements()
{
    Element onsite01(0, 1, [this](double kx, double y)
                     { return delta_SC; });

    Element onsite10(1, 0, [this](double kx, double y)
                     { return -delta_SC; });

    Delta_discrete_ky_onsite_elements = {onsite01, onsite10};
}

void ToyModelElementBased::set_Delta_Adjoint_discrete_ky_elements()
{
    Element onsite01(0, 1, [this](double kx, double y)
                     { return -delta_SC; });

    Element onsite10(1, 0, [this](double kx, double y)
                     { return delta_SC; });

    Delta_Adjoint_discrete_ky_onsite_elements = {onsite01, onsite10};
}

void ToyModelElementBased::set_mHmkT_discrete_ky_elements()
{
    Element onsite00(0, 0, [this](double kx, double y)
                     { return -Ek(kx) - HBz() + mu; });

    Element onsite01(0, 1, [this](double kx, double y)
                     { return -std::conj(HBx_iy()) + std::conj(HRSO(kx)); });

    Element onsite10(1, 0, [this](double kx, double y)
                     { return -HBx_iy() + HRSO(kx); });

    Element onsite11(1, 1, [this](double kx, double y)
                     { return -Ek(kx) + HBz() + mu; });

    mHmkT_discrete_ky_onsite_elements = {onsite00, onsite01, onsite10, onsite11};

    Element hopping_p00(0, 0, [this](double kx, double y)
                        { return ty; });

    Element hopping_p01(0, 1, [this](double kx, double y)
                        { return -1i * 0.5 * delta_RSO_y; });

    Element hopping_p10(1, 0, [this](double kx, double y)
                        { return -1i * 0.5 * delta_RSO_y; });

    Element hopping_p11(1, 1, [this](double kx, double y)
                        { return ty; });

    mHmkT_discrete_ky_hopping_p_elements = {hopping_p00, hopping_p01, hopping_p10, hopping_p11};

    Element hopping_m00(0, 0, [this](double kx, double y)
                        { return ty; });

    Element hopping_m01(0, 1, [this](double kx, double y)
                        { return 1i * 0.5 * delta_RSO_y; });

    Element hopping_m10(1, 0, [this](double kx, double y)
                        { return 1i * 0.5 * delta_RSO_y; });

    Element hopping_m11(1, 1, [this](double kx, double y)
                        { return ty; });

    mHmkT_discrete_ky_hopping_m_elements = {hopping_m00, hopping_m01, hopping_m10, hopping_m11};
}

void ToyModelElementBased::set_Hk_discrete_elements()
{
    Element onsite00(0, 0, [this](double x, double y)
                     { return Ek() + HBz() - mu; });

    Element onsite01(0, 1, [this](double x, double y)
                     { return HBx_iy(); });

    Element onsite10(1, 0, [this](double x, double y)
                     { return std::conj(HBx_iy()); });

    Element onsite11(1, 1, [this](double x, double y)
                     { return Ek() - HBz() - mu; });

    Hk_discrete_onsite_elements = {onsite00, onsite01, onsite10, onsite11};

    Element hopping_xp00(0, 0, [this](double x, double y)
                         { return -tx; });

    Element hopping_xp01(0, 1, [this](double x, double y)
                         { return 0.5 * delta_RSO_x; });

    Element hopping_xp10(1, 0, [this](double x, double y)
                         { return -0.5 * delta_RSO_x; });

    Element hopping_xp11(1, 1, [this](double x, double y)
                         { return -tx; });

    Hk_discrete_hopping_xp_elements = {hopping_xp00, hopping_xp01, hopping_xp10, hopping_xp11};

    Element hopping_xm00(0, 0, [this](double x, double y)
                         { return -tx; });

    Element hopping_xm01(0, 1, [this](double x, double y)
                         { return -0.5 * delta_RSO_x; });

    Element hopping_xm10(1, 0, [this](double x, double y)
                         { return 0.5 * delta_RSO_x; });

    Element hopping_xm11(1, 1, [this](double x, double y)
                         { return -tx; });

    Hk_discrete_hopping_xm_elements = {hopping_xm00, hopping_xm01, hopping_xm10, hopping_xm11};

    Element hopping_yp00(0, 0, [this](double x, double y)
                         { return -ty; });

    Element hopping_yp01(0, 1, [this](double x, double y)
                         { return -1i * 0.5 * delta_RSO_y; });

    Element hopping_yp10(1, 0, [this](double x, double y)
                         { return -1i * 0.5 * delta_RSO_y; });

    Element hopping_yp11(1, 1, [this](double x, double y)
                         { return -ty; });

    Hk_discrete_hopping_yp_elements = {hopping_yp00, hopping_yp01, hopping_yp10, hopping_yp11};

    Element hopping_ym00(0, 0, [this](double x, double y)
                         { return -ty; });

    Element hopping_ym01(0, 1, [this](double x, double y)
                         { return 1i * 0.5 * delta_RSO_y; });

    Element hopping_ym10(1, 0, [this](double x, double y)
                         { return 1i * 0.5 * delta_RSO_y; });

    Element hopping_ym11(1, 1, [this](double x, double y)
                         { return -ty; });

    Hk_discrete_hopping_ym_elements = {hopping_ym00, hopping_ym01, hopping_ym10, hopping_ym11};
}

void ToyModelElementBased::set_Delta_discrete_elements()
{
    Element onsite01(0, 1, [this](double x, double y)
                     { return delta_SC; });

    Element onsite10(1, 0, [this](double x, double y)
                     { return -delta_SC; });

    Delta_discrete_onsite_elements = {onsite01, onsite10};
}

void ToyModelElementBased::set_Delta_Adjoint_discrete_elements()
{
    Element onsite01(0, 1, [this](double x, double y)
                     { return -delta_SC; });

    Element onsite10(1, 0, [this](double x, double y)
                     { return delta_SC; });

    Delta_Adjoint_discrete_onsite_elements = {onsite01, onsite10};
}

void ToyModelElementBased::set_mHmkT_discrete_elements()
{
    Element onsite00(0, 0, [this](double x, double y)
                     { return -Ek() - HBz() + mu; });

    Element onsite01(0, 1, [this](double x, double y)
                     { return -std::conj(HBx_iy()); });

    Element onsite10(1, 0, [this](double x, double y)
                     { return -HBx_iy(); });

    Element onsite11(1, 1, [this](double x, double y)
                     { return -Ek() + HBz() + mu; });

    mHmkT_discrete_onsite_elements = {onsite00, onsite01, onsite10, onsite11};

    Element hopping_xp00(0, 0, [this](double x, double y)
                         { return tx; });

    Element hopping_xp01(0, 1, [this](double x, double y)
                         { return -0.5 * delta_RSO_x; });

    Element hopping_xp10(1, 0, [this](double x, double y)
                         { return 0.5 * delta_RSO_x; });

    Element hopping_xp11(1, 1, [this](double x, double y)
                         { return tx; });

    mHmkT_discrete_hopping_xp_elements = {hopping_xp00, hopping_xp01, hopping_xp10, hopping_xp11};

    Element hopping_xm00(0, 0, [this](double x, double y)
                         { return tx; });

    Element hopping_xm01(0, 1, [this](double x, double y)
                         { return 0.5 * delta_RSO_x; });

    Element hopping_xm10(1, 0, [this](double x, double y)
                         { return -0.5 * delta_RSO_x; });

    Element hopping_xm11(1, 1, [this](double x, double y)
                         { return tx; });

    mHmkT_discrete_hopping_xm_elements = {hopping_xm00, hopping_xm01, hopping_xm10, hopping_xm11};

    Element hopping_yp00(0, 0, [this](double x, double y)
                         { return ty; });

    Element hopping_yp01(0, 1, [this](double x, double y)
                         { return -1i * 0.5 * delta_RSO_y; });

    Element hopping_yp10(1, 0, [this](double x, double y)
                         { return -1i * 0.5 * delta_RSO_y; });

    Element hopping_yp11(1, 1, [this](double x, double y)
                         { return ty; });

    mHmkT_discrete_hopping_yp_elements = {hopping_yp00, hopping_yp01, hopping_yp10, hopping_yp11};

    Element hopping_ym00(0, 0, [this](double x, double y)
                         { return ty; });

    Element hopping_ym01(0, 1, [this](double x, double y)
                         { return 1i * 0.5 * delta_RSO_y; });

    Element hopping_ym10(1, 0, [this](double x, double y)
                         { return 1i * 0.5 * delta_RSO_y; });

    Element hopping_ym11(1, 1, [this](double x, double y)
                         { return ty; });

    mHmkT_discrete_hopping_ym_elements = {hopping_ym00, hopping_ym01, hopping_ym10, hopping_ym11};
}