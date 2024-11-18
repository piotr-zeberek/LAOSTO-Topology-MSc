#include "System2D.h"
#include "utils.h"

#include <iostream>
#include <fstream>
#include <thread>

#define ARMA_USE_ARPACK 1
#define ARMA_USE_SUPERLU 1

#include <armadillo>

System2D::System2D(const H2D &H, const Parameters &p) : _H(H), _p(p)
{
    _n_bands = _H({0, 0}, _p).rows();
}

Eigen::VectorXd System2D::eigenvals(const Point2D &k)
{
    _SAES.compute(_H(k, _p), Eigen::EigenvaluesOnly);
    return _SAES.eigenvalues();
};

Eigen::MatrixXcd System2D::eigenvecs(const Point2D &k)
{
    _SAES.compute(_H(k, _p));
    return _SAES.eigenvectors();
};

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> System2D::eigen(const Point2D &k)
{
    _SAES.compute(_H(k, _p));
    return std::make_tuple(_SAES.eigenvalues(), _SAES.eigenvectors());
};

void System2D::printBandStructure(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
    {
        for (auto ky : ky_vec)
        {
            _SAES.compute(_H({kx, ky}, _p), Eigen::EigenvaluesOnly);
            output_file << kx << " " << ky << " " << _SAES.eigenvalues().transpose() / meV2au(1.0) << "\n";
        }
    }
}

void System2D::printBandStructureSlice(const std::string &output_filename, const Eigen::VectorXd &k_vec, int axis, double k0)
{
    std::ofstream output_file(output_filename);

    if (axis == 0)
    {
        for (auto kx : k_vec)
        {
            _SAES.compute(_H({kx, k0}, _p), Eigen::EigenvaluesOnly);
            output_file << kx << " " << _SAES.eigenvalues().transpose() / meV2au(1.0) << "\n";
        }
    }
    else if (axis == 1)
    {
        for (auto ky : k_vec)
        {
            _SAES.compute(_H({k0, ky}, _p), Eigen::EigenvaluesOnly);
            output_file << ky << " " << _SAES.eigenvalues().transpose() / meV2au(1.0) << "\n";
        }
    }
    else
    {
        std::cerr << "Invalid axis. Choose 0 or 1." << std::endl;
    }
}

void System2D::printBandStructureDiscreteky(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky, bool use_arma)
{
    std::ofstream output_file(output_filename);

    Eigen::VectorXd evals(_n_bands * n_ky);

    for (auto kx : kx_vec)
    {
        if (use_arma)
        {
            evals = H_discrete_ky_eigendecomposiiton(kx, n_ky).first;
            std::sort(evals.data(), evals.data() + evals.size());
        }
        else
        {
            evals = _SAES.compute(H_discrete_ky(kx, n_ky), Eigen::EigenvaluesOnly).eigenvalues();
        }

        output_file << kx << " " << evals.transpose() / meV2au(1.0) << std::endl;
    }
}

void System2D::printProbDenDiscrete(const std::string &output_filename, std::size_t nk_x, std::size_t nk_y, double E)
{
    std::ofstream output_file(output_filename);

    auto [evals, evecs] = H_discrete_eigendecomposiiton(nk_x, nk_y, false, 1, E);
    Eigen::VectorXd prob_den_orbitals = evecs.col(0).cwiseAbs2();

    output_file << "# E = " << E / meV2au(1) << " eval = " << evals(0) / meV2au(1) << std::endl;

    for (auto i = 0; i < nk_x; ++i)
    {
        for (auto j = 0; j < nk_y; ++j)
        {
            output_file << i << " " << j << " " << prob_den_orbitals.segment((i * nk_y + j) * _n_bands, _n_bands).sum() << std::endl;
        }
    }
}

void System2D::printOrdinaryGap(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
    {
        for (auto ky : ky_vec)
        {
            _SAES.compute(_H({kx, ky}, _p), Eigen::EigenvaluesOnly);
            output_file << kx << " " << ky << " " << (_SAES.eigenvalues().tail(_n_bands / 2) - _SAES.eigenvalues().head(_n_bands / 2).reverse()).transpose() / meV2au(1.0) << "\n";
        }
    }
}

void System2D::printOrdinaryGapAlongContour(const std::string &output_filename, const std::vector<Point2D> &contour)
{
    std::ofstream output_file(output_filename);

    for (auto k : contour)
    {
        _SAES.compute(_H(k, _p), Eigen::EigenvaluesOnly);
        output_file << k.x() << " " << k.y() << " " << (_SAES.eigenvalues().tail(_n_bands / 2) - _SAES.eigenvalues().head(_n_bands / 2).reverse()).transpose() / meV2au(1.0) << "\n";
    }
}

void System2D::printAbsDeltaAlongContour(const std::string &output_filename, const std::vector<Point2D> &contour)
{
    std::ofstream output_file(output_filename);

    for (auto k : contour)
    {
        output_file << k.x() << " " << k.y() << " " << calcAbsDelta(k).transpose() / meV2au(1) / meV2au(1) << "\n";
    }
}

// void System2D::printDeltaFromUnitaryTransformation(const std::string &delta_filename, const std::string &DT_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec)
// {
//     std::ofstream delta_file(delta_filename);
//     std::ofstream DT_file(DT_filename);

//     std::size_t half_bands = _n_bands / 2;

//     Eigen::MatrixXcd H(_n_bands, _n_bands);
//     Eigen::MatrixXcd U_tl(half_bands, half_bands);
//     Eigen::MatrixXcd U_br(half_bands, half_bands);
//     Eigen::MatrixXcd U = Eigen::MatrixXcd::Zero(_n_bands, _n_bands);

//     Eigen::VectorXd evals(half_bands);

//     // std::vector<std::size_t> indices(half_bands, 0);
//     // std::vector<std::size_t> inv_indices(half_bands, 0);

//     // for (auto kx : kx_vec)
//     // {
//     //     for (auto ky : ky_vec)
//     //     {
//     //         H = _H({kx, ky}, _p);

//     //         auto U_tl = _SAES.compute(H.topLeftCorner(half_bands, half_bands)).eigenvectors();
//     //         auto U_br = _SAES.compute(H.bottomRightCorner(half_bands, half_bands)).eigenvectors().rowwise().reverse();
//     //         evals = -_SAES.eigenvalues();

//     //         _SAES.compute(H, Eigen::EigenvaluesOnly);

//     //         // reset idx
//     //         for (auto i = 0; i < half_bands; ++i)
//     //         {
//     //             indices[i] = i;
//     //         }

//     //         std::sort(indices.begin(), indices.end(), [evals](std::size_t i1, std::size_t i2)
//     //                   { return std::abs(evals(i1)) <= std::abs(evals(i2)); });

//     //         for (auto i = 0; i < half_bands; ++i)
//     //         {
//     //             inv_indices[indices[i]] = i;
//     //         }

//     //         // for (auto i = 0; i < half_bands; ++i)
//     //         // {
//     //         //     std::size_t idx = 0;
//     //         //     U_tl.col(i).cwiseAbs2().maxCoeff(&idx);
//     //         //     U_tl.col(i) /= U_tl.col(i)(idx);
//     //         //     U_tl.col(i).normalize();

//     //         //     U_br.col(i).cwiseAbs2().maxCoeff(&idx);
//     //         //     U_br.col(i) /= U_br.col(i)(idx);
//     //         //     U_br.col(i).normalize();
//     //         // }

//     //         // std::cout << (U.adjoint() * sys._H({0, 0.0}, p) * U)/meV2au(1) << std::endl;

//     //         U.topLeftCorner(half_bands, half_bands) = U_tl;
//     //         U.bottomRightCorner(half_bands, half_bands) = U_br;

//     //         for (auto i = 0; i < half_bands; ++i)
//     //         {
//     //             abs2_deltas(i) = _SAES.eigenvalues()(half_bands + inv_indices[i]) * _SAES.eigenvalues()(half_bands + inv_indices[i]) - evals(i) * evals(i);
//     //         }

//     //         // abs2_deltas += _SAES.eigenvalues().tail(half_bands).cwiseAbs2();

//     //         output_file << kx << " " << ky << " " << abs2_deltas.transpose() / meV2au(1) / meV2au(1) << "\n";
//     //     }
//     // }

//     delta_file << "# kx ky Re(d0_00) Im(d0_00) Re(d0_01) Im(d0_01) Re(d0_10) Im(d0_10) Re(d0_11) Im(d0_11) ... \n";
//     DT_file << "# kx ky Re(det(d0)) Im(det(d0)) Re(tr(d0)) Im(tr(d0)) ... \n";

//     double eps = 1e-3;

//     auto setPhase = [half_bands](Eigen::MatrixXcd &U_tl, Eigen::MatrixXcd &U_br)
//     {
//         std::size_t idx = 0;

//         for (auto i = 0; i < half_bands; ++i)
//         {
//             (U_tl.col(i).cwiseAbs() + U_br.col(i).cwiseAbs()).maxCoeff(&idx);
//             U_tl.col(i) /= U_tl.col(i)(idx);
//             U_br.col(i) /= U_br.col(i)(idx);

//             U_tl.col(i).normalize();
//             U_br.col(i).normalize();
//         }
//     };

//     std::ofstream test_file("data/test.dat");
//     std::ofstream test2_file("data/test2.dat");
//     std::ofstream test3_file("data/test3.dat");
//     std::ofstream test4_file("data/test4.dat");

//     Eigen::VectorXcd ref = Eigen::VectorXcd::Constant(half_bands, Eigen::dcomplex(std::sqrt(1.0 / half_bands), 0));
//     Eigen::dcomplex inner_product{};

//     for (auto kx : kx_vec)
//     {
//         for (auto ky : ky_vec)
//         {
//             H = _H({kx, ky}, _p);

//             test_file << kx << " " << ky << " ";
//             Eigen::MatrixXcd U_tl = _SAES.compute(H.topLeftCorner(half_bands, half_bands)).eigenvectors();
//             test_file << _SAES.eigenvalues().transpose() / meV2au(1) << " ";
//             Eigen::MatrixXcd U_br = _SAES.compute(H.bottomRightCorner(half_bands, half_bands)).eigenvectors().rowwise().reverse();
//             test_file << _SAES.eigenvalues().transpose().reverse() / meV2au(1) << std::endl;

//             // // idx th element
//             // std::size_t idx = 5;
//             // for (auto i = 0; i < half_bands; ++i)
//             // {
//             //     U_tl.col(i) /= U_tl.col(i)(idx);
//             //     U_br.col(i) /= U_br.col(i)(idx);

//             //     U_tl.col(i).normalize();
//             //     U_br.col(i).normalize();
//             // }

//             // // max abs sum
//             // std::size_t idx = 0;
//             // for (auto i = 0; i < half_bands; ++i)
//             // {
//             //     (U_tl.col(i).cwiseAbs2() + U_br.col(i).cwiseAbs2()).maxCoeff(&idx);
//             //     U_tl.col(i) /= U_tl.col(i)(idx);
//             //     U_br.col(i) /= U_br.col(i)(idx);

//             //     U_tl.col(i).normalize();
//             //     U_br.col(i).normalize();
//             // }

//             // // max abs of outer product
//             // std::size_t idx_tl = 0;
//             // std::size_t idx_br = 0;
//             // for (auto i = 0; i < half_bands; ++i)
//             // {
//             //     (U_tl.col(i) * U_br.col(i).transpose()).cwiseAbs2().maxCoeff(&idx_tl, &idx_br);

//             //     U_tl.col(i) /= U_tl.col(i)(idx_tl);
//             //     U_br.col(i) /= U_br.col(i)(idx_br);

//             //     U_tl.col(i).normalize();
//             //     U_br.col(i).normalize();
//             // }

//             // // phase of the inner product
//             // std::size_t idx = 0;
//             // for (auto i = 0; i < half_bands; ++i)
//             // {
//             //     double phase = std::arg(U_tl.col(i).dot(U_br.col(i)));

//             //     U_tl.col(i) *= std::exp(-1i * phase/2.0);
//             //     U_br.col(i) *= std::exp(-1i * phase/2.0);
//             // }

//             // // phase of the all inner products
//             // double phase = (U_tl.adjoint()*U_br).cwiseArg().sum();
//             // U_tl *= std::exp(-1i * phase/2.0);1
//             // U_br *= std::exp(-1i * phase/2.0);

//             // // first vector of U_tl as reference
//             // Eigen::VectorXcd ref = U_tl.col(0);
//             // Eigen::dcomplex inner_product;

//             // inner_product = ref.dot(U_br.col(0));
//             // U_br.col(0) *= inner_product / std::abs(inner_product);

//             // for (auto i = 1; i < half_bands; ++i)
//             // {
//             //     inner_product = ref.dot(U_tl.col(i));
//             //     U_tl.col(i) *= inner_product / std::abs(inner_product);

//             //     inner_product = ref.dot(U_br.col(i));
//             //     U_br.col(i) *= inner_product / std::abs(inner_product);
//             // }

//             // custom reference
//             for (auto i = 0; i < half_bands; ++i)
//             {
//                 U_tl.col(i) *= U_tl.col(i).dot(ref);
//                 U_br.col(i) *= U_br.col(i).dot(ref);

//                 U_tl.col(i).normalize();
//                 U_br.col(i).normalize();
//             }

//             Eigen::MatrixXcd delta_matrix = U_tl.adjoint() * H.topRightCorner(half_bands, half_bands) * U_br;
//             Eigen::VectorXcd delta_matrix_vector_view = Eigen::Map<Eigen::VectorXcd>(delta_matrix.data(), delta_matrix.size());

//             test2_file << kx << " " << ky << " " << delta_matrix_vector_view.cwiseAbs().transpose() / meV2au(1) << std::endl;
//             test3_file << kx << " " << ky << " ";
//             test4_file << kx << " " << ky << " " << std::sqrt(0.5 * (delta_matrix.adjoint() * delta_matrix).trace().real()) / meV2au(1) << " "
//                        << std::sqrt(0.5 * (delta_matrix.adjoint() * delta_matrix).trace().imag()) / meV2au(1) << std::endl;

//             delta_file << kx << " " << ky << " ";
//             DT_file << kx << " " << ky << " ";

//             for (auto i = 0; i < half_bands; i = i + 2)
//             {
//                 Eigen::Matrix2cd d = delta_matrix.block(i, i, 2, 2);

//                 for (auto j = 0; j < 4; ++j)
//                 {
//                     delta_file << d(j / 2, j % 2).real() / meV2au(1) << " " << d(j / 2, j % 2).imag() / meV2au(1) << " ";
//                 }

//                 DT_file << d.determinant().real() / meV2au(1) / meV2au(1) << " " << d.determinant().imag() / meV2au(1) / meV2au(1) << " "
//                         << d.trace().real() / meV2au(1) << " " << d.trace().imag() / meV2au(1) << " ";

//                 test3_file << (d.adjoint() * d).trace().real() / meV2au(1) / meV2au(1) << " " << (d.adjoint() * d).trace().imag() / meV2au(1) / meV2au(1) << " ";
//             }

//             delta_file << "\n";
//             DT_file << "\n";
//             test3_file << "\n";
//         }
//     }
// }

void System2D::printDeltaFromUnitaryTransformation(const std::string &delta_filename, const std::string &DT_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec)
{
    std::ofstream delta_file(delta_filename);
    std::ofstream DT_file(DT_filename);

    std::size_t half_bands = _n_bands / 2;
    std::size_t orbital_bands = half_bands / 2;

    Eigen::MatrixXcd H(_n_bands, _n_bands);
    Eigen::MatrixXcd Io = Eigen::MatrixXcd::Identity(orbital_bands, orbital_bands);

    // chiral symmetry U
    //  Eigen::MatrixXcd U = (kron(sz, kron(Io, s0)) + kron(sx, kron(Io, sy))) / std::sqrt(2.0);

    // PHS U
    Eigen::MatrixXcd C = kron(sx, kron(Io, s0));
    Eigen::MatrixXcd D = Eigen::MatrixXcd::Zero(_n_bands, _n_bands);
    _SAES.compute(C);

    for (auto i = 0; i < _n_bands; ++i)
    {
        D(i, i) = _SAES.eigenvalues()(i) > 0 ? 1 / std::sqrt(_SAES.eigenvalues()(i)) : Eigen::dcomplex(0, 1 / std::sqrt(-_SAES.eigenvalues()(i)));
    }
    Eigen::MatrixXcd U = _SAES.eigenvectors() * D;

    for (auto kx : kx_vec)
    {
        for (auto ky : ky_vec)
        {
            H = _H({kx, ky}, _p);

            H = (U.adjoint() * H * U).eval();

            // std::cout << H / meV2au(1) << std::endl;

            Eigen::MatrixXcd delta_matrix = H.topRightCorner(half_bands, half_bands);
            Eigen::VectorXcd delta_matrix_vector_view = Eigen::Map<Eigen::VectorXcd>(delta_matrix.data(), delta_matrix.size());

            delta_file << kx << " " << ky << " " << delta_matrix_vector_view.cwiseAbs().transpose() / meV2au(1) << std::endl;

            Eigen::dcomplex det = delta_matrix.determinant();
            Eigen::dcomplex tr = delta_matrix.trace();

            DT_file << kx << " " << ky << " " << det.real() / meV2au(1) / meV2au(1) << " " << det.imag() / meV2au(1) / meV2au(1) << " "
                    << tr.real() / meV2au(1) << " " << tr.imag() / meV2au(1) << std::endl;
        }
    }
}

void System2D::printAbsDelta(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
    {
        for (auto ky : ky_vec)
        {
            output_file << kx << " " << ky << " " << calcAbsDelta({kx, ky}).transpose() / meV2au(1) / meV2au(1) << std::endl;
        }
    }
}

Eigen::VectorXd System2D::calcAbsDelta(const Point2D &k)
{
    std::size_t half_bands = _n_bands / 2;

    Eigen::VectorXd abs_deltas(_n_bands);

    Eigen::MatrixXcd H = _H(k, _p);

    // Eigen::VectorXd Ep = _SAES.compute(H.topLeftCorner(half_bands, half_bands), Eigen::EigenvaluesOnly).eigenvalues();
    // Eigen::VectorXd Em = _SAES.compute(H.bottomRightCorner(half_bands, half_bands), Eigen::EigenvaluesOnly).eigenvalues();

    // _SAES.compute(H, Eigen::EigenvaluesOnly);

    // Eigen::VectorXd diffs(_n_bands);
    // std::size_t idx = 0;

    // for (auto i = 0; i < half_bands; ++i)
    // {
    //     (_SAES.eigenvalues().cwiseAbs2().array() - Ep(i) * Ep(i)).cwiseAbs().minCoeff(&idx);
    //     abs_deltas(i) = _SAES.eigenvalues()(idx) * _SAES.eigenvalues()(idx) - Ep(i) * Ep(i);
    // }

    // for (auto i = 0; i < half_bands; ++i)
    // {
    //     (_SAES.eigenvalues().cwiseAbs2().array() - Em(i) * Em(i)).cwiseAbs().minCoeff(&idx);
    //     abs_deltas(half_bands + i) = _SAES.eigenvalues()(idx) * _SAES.eigenvalues()(idx) - Em(i) * Em(i);
    // }

    // return abs_deltas;

    // finding zeros of determinants
    {
        std::size_t half_bands = _n_bands / 2;

        // result vector
        Eigen::VectorXd abs_deltas(_n_bands);

        Eigen::MatrixXcd H = _H(k, _p);

        Eigen::MatrixXcd U_tl = _SAES.compute(H.topLeftCorner(half_bands, half_bands)).eigenvectors();
        Eigen::MatrixXd Ep = _SAES.eigenvalues().asDiagonal();

        Eigen::MatrixXcd U_br = _SAES.compute(H.bottomRightCorner(half_bands, half_bands)).eigenvectors();
        Eigen::MatrixXd Em = _SAES.eigenvalues().asDiagonal();

        Eigen::MatrixXcd U = Eigen::MatrixXcd::Zero(_n_bands, _n_bands);
        U.topLeftCorner(half_bands, half_bands) = U_tl;
        U.bottomRightCorner(half_bands, half_bands) = U_br;

        H = (U.adjoint() * H * U).eval();

        const Eigen::MatrixXcd &D = H.topRightCorner(half_bands, half_bands);
        const Eigen::MatrixXcd &D_hc = H.bottomLeftCorner(half_bands, half_bands);

        Eigen::MatrixXcd I = Eigen::MatrixXcd::Identity(half_bands, half_bands);

        auto det_Ep_inv = [Ep, Em, D, D_hc, I](long double lambda) -> long double
        {
            return (Em - lambda * I - D_hc * (Ep - lambda * I).inverse() * D).determinant().real();
        };

        auto det_Em_inv = [Ep, Em, D, D_hc, I](long double lambda) -> long double
        {
            return (Ep - lambda * I - D * (Em - lambda * I).inverse() * D_hc).determinant().real();
        };

        auto check_for_changing_range = [](long double detl, long double detr) -> int
        {
            if (detl * detr > 0)
            {
                if (detl + detr > 0)
                {
                    if (detr - detl > 0)
                    {
                        return -1;
                    }
                    else
                    {
                        return 1;
                    }
                }
                else
                {
                    if (detr - detl > 0)
                    {
                        return 1;
                    }
                    else
                    {
                        return -1;
                    }
                }
            }

            return 0;
        };

        long double lam{};
        long double laml{};
        long double lamr{};

        long double det{};
        long double detl{};
        long double detr{};

        long double dlam = meV2au(1e-6);
        long double eps = meV2au(1e-12);

        const std::size_t max_iter = 100;

        for (auto i = 0; i < half_bands; ++i)
        {
            lam = Ep(i, i);

            laml = lam - dlam;
            lamr = lam + dlam;

            detl = det_Em_inv(laml);
            detr = det_Em_inv(lamr);

            int change_range = check_for_changing_range(detl, detr);

            switch (change_range)
            {
            case -1:
                do
                {
                    laml -= dlam;
                    detl = det_Em_inv(laml);
                } while (detl * detr > 0);

                lamr = lam + dlam;
                detr = det_Em_inv(lamr);
                break;

            case 1:
                do
                {
                    lamr += dlam;
                    detr = det_Em_inv(lamr);
                } while (detl * detr > 0);

                laml = lam - dlam;
                detl = det_Em_inv(laml);
                break;

            default:
                break;
            }

            for (std::size_t iter = 0; iter < max_iter && std::abs(lamr - laml) > eps; ++iter)
            {

                lam = (detl * lamr - detr * laml) / (detl - detr);
                det = det_Em_inv(lam);

                if (detl * det < 0)
                {
                    lamr = lam;
                    detr = det;
                }
                else
                {
                    laml = lam;
                    detl = det;
                }
            }

            abs_deltas(i) = lam * lam - Ep(i, i) * Ep(i, i);
        }

        for (auto i = 0; i < half_bands; ++i)
        {
            lam = Em(i, i);

            laml = lam - dlam;
            lamr = lam + dlam;

            detl = det_Ep_inv(laml);
            detr = det_Ep_inv(lamr);

            int change_range = check_for_changing_range(detl, detr);

            switch (change_range)
            {
            case -1:
                do
                {
                    laml -= dlam;
                    detl = det_Ep_inv(laml);
                } while (detl * detr > 0);

                lamr = lam + dlam;
                detr = det_Ep_inv(lamr);
                break;

            case 1:
                do
                {
                    lamr += dlam;
                    detr = det_Ep_inv(lamr);
                } while (detl * detr > 0);

                laml = lam - dlam;
                detl = det_Ep_inv(laml);
                break;

            default:
                break;
            }

            for (std::size_t iter = 0; iter < max_iter && std::abs(lamr - laml) > eps; ++iter)
            {
                lam = (detl * lamr - detr * laml) / (detl - detr);
                det = det_Ep_inv(lam);

                if (detl * det < 0)
                {
                    lamr = lam;
                    detr = det;
                }
                else
                {
                    laml = lam;
                    detl = det;
                }
            }

            abs_deltas(i + half_bands) = lam * lam - Em(i, i) * Em(i, i);
        }

        return abs_deltas;
    }
}

std::vector<std::vector<Point2D>> System2D::findFSContours(double E, double dk, double eps, const Point2D kx_range, std::size_t n_kx)
{
    std::vector<std::vector<Point2D>> FS_contours;

    // search for kx such that E is an eigenvalue of H(kx, 0)
    Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_kx, kx_range(0), kx_range(1));

    _SAES.compute(_H({kx_vec(0), 0.0}, _p), Eigen::EigenvaluesOnly);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> SAES_next;

    std::vector<std::size_t> kx_idx; // indices of kx_vec where the contour is found
    std::vector<std::size_t> E_idx;  // indices of bands creating the contour

    for (auto ie = 0; ie < _n_bands; ++ie) // contour index
    {
        for (auto i = (kx_idx.empty() ? 0 : kx_idx.back()); i < kx_vec.size() - 1; ++i)
        {
            switch (i % 2)
            {
            case 0:
                SAES_next.compute(_H({kx_vec(i + 1), 0.0}, _p), Eigen::EigenvaluesOnly);
                break;
            case 1:
                _SAES.compute(_H({kx_vec(i + 1), 0.0}, _p), Eigen::EigenvaluesOnly);
                break;
            }

            if ((_SAES.eigenvalues()(ie) - E) * (SAES_next.eigenvalues()(ie) - E) < 0)
            {
                // std::cout << "contour " << ic << " found at kx = " << kx_vec(i) << std::endl;
                kx_idx.push_back(i);
                E_idx.push_back(ie);
                break;
            }
        }
    }

    auto R_mat = [](double phi) -> Eigen::Matrix2d
    {
        Eigen::Matrix2d R;
        R << std::cos(phi), -std::sin(phi), std::sin(phi), std::cos(phi);
        return R;
    };

    // begin looking for contours in -ky direction
    for (auto ic = 0; ic < kx_idx.size(); ++ic)
    {
        // find exact start kx with regula falsi
        double kxl = kx_vec[kx_idx[ic]];
        double kxr = kx_vec[kx_idx[ic] + 1];

        double El = _SAES.compute(_H({kxl, 0.0}, _p), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);
        double Er = _SAES.compute(_H({kxr, 0.0}, _p), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);

        double kx_exact = 0.0;
        double E_exact = 0.0;

        do
        {
            kx_exact = (El * kxr - Er * kxl) / (El - Er);
            E_exact = _SAES.compute(_H({kx_exact, 0.0}, _p), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);

            // if (std::abs(E_exact - E) < eps)
            //     break;

            if (El * E_exact < 0)
            {
                kxr = kx_exact;
                Er = E_exact;
            }
            else
            {
                kxl = kx_exact;
                El = E_exact;
            }

        } while (std::abs(E_exact - E) > eps);

        // std::cout << "contour " << ic << " found at kx_exact = " << kx_exact << std::endl;

        // first point of the contour
        std::vector<Point2D> contour;
        contour.push_back(Point2D{kx_exact, 0.0});
        contour.reserve(1e5); // try to avoid reallocations

        // find next points of the contour
        double phil = 0.0;
        double phir = 0.0;

        Eigen::Matrix2d Rl;
        Eigen::Matrix2d Rr;

        Point2D k_diff;
        Point2D kl;
        Point2D kr;

        double phi_next = 0.0;
        Eigen::Matrix2d R_next;
        Point2D k_next = {0.0, 0.0};
        double E_next = 0.0;

        do
        {
            phil = M_PI_4;
            phir = -M_PI_4;

            Rl = R_mat(phil);
            Rr = Rl;
            Rr(0, 1) *= -1;
            Rr(1, 0) *= -1;

            k_diff = contour.size() > 1 ? contour.end()[-1] - contour.end()[-2] : Point2D{0.0, -dk};

            kl = contour.end()[-1] + Rl * k_diff;
            kr = contour.end()[-1] + Rr * k_diff;

            El = _SAES.compute(_H(kl, _p), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);
            Er = _SAES.compute(_H(kr, _p), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);

            // find next k - regula falsi with phi
            do
            {
                phi_next = (El * phir - Er * phil) / (El - Er);
                R_next = R_mat(phi_next);
                k_next = contour.end()[-1] + R_next * k_diff;
                E_next = _SAES.compute(_H(k_next, _p), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);

                // if (std::abs(E_next - E) < eps)
                //     break;

                if (El * E_exact < 0)
                {
                    phir = phi_next;
                    Er = E_next;
                }
                else
                {
                    phil = phi_next;
                    El = E_next;
                }

            } while (std::abs(E_next - E) > eps);

            contour.push_back(k_next);
        } while (contour.size() < 3 || (contour.end()[-1] - contour[0]).norm() > dk);

        // close the contour
        if (contour.back().y() < 0)
            contour.back() = contour[0];
        else
            contour.push_back(contour[0]);

        FS_contours.push_back(contour);
    }

    return FS_contours;
}

void System2D::printBerryCurvature(const std::string &output_filename, const H2D &dvxH, const H2D &dvyH, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec)
{
    std::ofstream output_file(output_filename);

    Eigen::VectorXd BC(_n_bands); // Berry curvature

    Eigen::MatrixXcd dvxH_mat(_n_bands, _n_bands);
    Eigen::MatrixXcd dvyH_mat(_n_bands, _n_bands);

    Eigen::dcomplex v1x, v1y, v2x, v2y;

    for (auto kx : kx_vec)
    {
        for (auto ky : ky_vec)
        {
            _SAES.compute(_H({kx, ky}, _p));

            dvxH_mat = dvxH({kx, ky}, _p);
            dvyH_mat = dvyH({kx, ky}, _p);

            BC.setZero();

            for (auto bi = 0; bi < _n_bands; ++bi) // bi - band index
            {

                for (auto bj = 0; bj < _n_bands; ++bj)
                {
                    if (bj == bi)
                        continue;

                    v1x = _SAES.eigenvectors().col(bi).dot(dvxH_mat * _SAES.eigenvectors().col(bj));
                    v1y = _SAES.eigenvectors().col(bi).dot(dvyH_mat * _SAES.eigenvectors().col(bj));

                    v2x = _SAES.eigenvectors().col(bj).dot(dvxH_mat * _SAES.eigenvectors().col(bi));
                    v2y = _SAES.eigenvectors().col(bj).dot(dvyH_mat * _SAES.eigenvectors().col(bi));

                    BC(bi) -= std::imag(v1x * v2y - v2x * v1y) / (_SAES.eigenvalues()(bi) - _SAES.eigenvalues()(bj)) / (_SAES.eigenvalues()(bi) - _SAES.eigenvalues()(bj));
                }
            }

            output_file << kx << " " << ky << " " << BC.transpose() * 1e6 << "\n";
        }
    }
}

// Multithreaded version - second try
Eigen::VectorXd System2D::calcChernNumbers(const Eigen::Vector2<std::size_t> &n, double kmax)
{
    if (n.x() % 2 != 0 || n.y() % 2 != 0)
    {
        std::cerr << "Number of k-points must be even in each direction. Returning empty vector." << std::endl;
        return Eigen::VectorXd();
    }

    // BZ mesh
    Eigen::VectorXd kx;
    Eigen::VectorXd ky;

    if (kmax >= M_PI)
    {
        kx = Eigen::VectorXd::LinSpaced(n.x(), -M_PI, M_PI);
        ky = Eigen::VectorXd::LinSpaced(n.y(), -M_PI, M_PI);

        // account for periodic boundary conditions
        kx(n.x() - 1) = kx(0);
        ky(n.y() - 1) = ky(0);
    }
    else
    {
        kx = Eigen::VectorXd::LinSpaced(n.x(), -kmax, kmax);
        ky = Eigen::VectorXd::LinSpaced(n.y(), -kmax, kmax);
    }

    return calcChernNumbersWithCustomGrid(kx, ky);
}

Eigen::VectorXd System2D::calcChernNumbersDenserCenter(std::size_t n_dense, std::size_t n_sparse, double k_val)
{
    // dense region
    Eigen::VectorXd kc = Eigen::VectorXd::LinSpaced(n_dense, -k_val, k_val);

    // without the dense region
    Eigen::VectorXd ko = Eigen::VectorXd::LinSpaced(n_sparse, -M_PI, M_PI);
    ko(ko.size() - 1) = -M_PI;

    // before and after the dense region
    std::size_t n = ((1.0 - k_val / M_PI) * n_sparse) / 2.0;
    if (n % 2 != 0)
        n++;

    Eigen::VectorXd k_before = Eigen::VectorXd::LinSpaced(n, -M_PI, -k_val);
    Eigen::VectorXd k_after = Eigen::VectorXd::LinSpaced(n, k_val, M_PI);
    k_after(k_after.size() - 1) = -M_PI;

    // kx next to the dense region
    std::size_t m = k_val / M_PI * n_sparse;
    if (m % 2 != 0)
        m++;
    Eigen::VectorXd kxn = Eigen::VectorXd::LinSpaced(m, -k_val, k_val);

    // std::cout << "n = " << n << ", m = " << m << std::endl;
    // std::cout << "dense region: " << kc(0) << "," << kc(0) << " - " << kc(kc.size() - 1) << "," << kc(kc.size() - 1) << std::endl;
    // std::cout << "top bar: " << k_before(0) << "," << ko(0) << " - " << k_before(k_before.size() - 1) << "," << ko(ko.size() - 1) << std::endl;
    // std::cout << "bottom bar: " << k_after(0) << "," << ko(0) << " - " << k_after(k_after.size() - 1) << "," << ko(ko.size() - 1) << std::endl;
    // std::cout << "left smaller bar: " << kxn(0) << "," << k_before(0) << " - " << kxn(kxn.size() - 1) << "," << k_before(k_before.size() - 1) << std::endl;
    // std::cout << "right smaller bar: " << kxn(0) << "," << k_after(0) << " - " << kxn(kxn.size() - 1) << "," << k_after(k_after.size() - 1) << std::endl;

    return calcChernNumbersWithCustomGrid(kc, kc)          // dense region
           + calcChernNumbersWithCustomGrid(k_before, ko)  // top bar
           + calcChernNumbersWithCustomGrid(k_after, ko)   // botom bar
           + calcChernNumbersWithCustomGrid(kxn, k_before) // left smaller bar
           + calcChernNumbersWithCustomGrid(kxn, k_after); // right smaller bar
}

Eigen::VectorXd System2D::calcChernNumbersWithCustomGrid(const Eigen::VectorXd &kx, const Eigen::VectorXd &ky)
{
    std::size_t nx = kx.size();
    std::size_t ny = ky.size();

    if (nx % 2 != 0 || ny % 2 != 0)
    {
        std::cerr << "Number of k-points must be even in each direction. Returning empty vector." << std::endl;
        return Eigen::VectorXd();
    }

    // BZ mesh
    Eigen::VectorXd kx_t = kx.head(nx / 2 + 1);
    Eigen::VectorXd kx_b = kx.tail(nx / 2);

    Eigen::VectorXd ky_l = ky.head(ny / 2 + 1);
    Eigen::VectorXd ky_r = ky.tail(ny / 2);

    // with rows storing and jumping in rows
    auto processGridSection = [this](const Eigen::VectorXd &kx, const Eigen::VectorXd &ky, Eigen::VectorXd &CN) -> void
    {
        CN.setZero();

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> SAES;

        std::size_t rows = kx.size();
        std::size_t cols = ky.size();

        Eigen::VectorX<Eigen::MatrixXcd> i_evecs(cols);  // Current row eigenvectors
        Eigen::VectorX<Eigen::MatrixXcd> ni_evecs(cols); // Next row eigenvectors

        Eigen::VectorX<Eigen::MatrixXcd> *i_evecs_ptr = &i_evecs;
        Eigen::VectorX<Eigen::MatrixXcd> *ni_evecs_ptr = &ni_evecs;

        // Compute eigenvectors for the first row
        for (auto j = 0; j < cols; ++j)
        {
            i_evecs(j) = SAES.compute(_H({kx(0), ky(j)}, _p)).eigenvectors();
        }

        Eigen::VectorXd BCdkxdky(_n_bands); // Berry curvature * dkx * dky
        BCdkxdky.setZero();

        Eigen::dcomplex yp, xp, ym, xm; // eigenvectors dot products, yp in y+ direction, xm in x- direction etc.

        for (auto i = 0; i < rows - 1; ++i)
        {
            (*ni_evecs_ptr)(0) = SAES.compute(_H({kx(i + 1), ky(0)}, _p)).eigenvectors();

            for (auto j = 0; j < cols - 1; ++j)
            {
                (*ni_evecs_ptr)(j + 1) = SAES.compute(_H({kx(i + 1), ky(j + 1)}, _p)).eigenvectors(); // contains eigenvectors and eigenvalues for next i, next j

                for (auto bi = 0; bi < _n_bands; ++bi) // bi - band index
                {
                    yp = (*i_evecs_ptr)(j).col(bi).dot((*i_evecs_ptr)(j + 1).col(bi));      // psi_{i,j} dot psi_{i,j+1}
                    xp = (*i_evecs_ptr)(j + 1).col(bi).dot((*ni_evecs_ptr)(j + 1).col(bi)); // psi_{i,j+1} dot psi_{i+1,j+1}
                    ym = (*ni_evecs_ptr)(j + 1).col(bi).dot((*ni_evecs_ptr)(j).col(bi));    // psi_{i+1,j+1} dot psi_{i+1,j}
                    xm = (*ni_evecs_ptr)(j).col(bi).dot((*i_evecs_ptr)(j).col(bi));         // psi_{i+1,j} dot psi_{i,j}

                    BCdkxdky(bi) = -std::arg(yp * xp * ym * xm); // std::arg uses std::atan2, maybe negative?
                    // BCdkxdky(bi) = -std::log(yp * xp * ym * xm).imag(); // std::arg uses std::atan2, maybe negative?
                }

                CN += BCdkxdky;
            }

            std::swap(i_evecs_ptr, ni_evecs_ptr);
        }

        CN /= (2.0 * M_PI);
    };

    Eigen::VectorXd CN_tl(_n_bands);
    Eigen::VectorXd CN_tr(_n_bands);
    Eigen::VectorXd CN_bl(_n_bands);
    Eigen::VectorXd CN_br(_n_bands);

    std::thread tl(processGridSection, kx_t, ky_l, std::ref(CN_tl));
    std::thread tr(processGridSection, kx_t, ky_r, std::ref(CN_tr));
    std::thread bl(processGridSection, kx_b, ky_l, std::ref(CN_bl));
    // std::thread br(processGridSection, kx_b, ky_r, std::ref(CN_br));

    // bottom right using main thread
    processGridSection(kx_b, ky_r, CN_br);

    tl.join();
    tr.join();
    bl.join();
    // br.join();

    return CN_tl + CN_tr + CN_bl + CN_br;
}

// w zasadzie dziala, ale jest wolne i dzwoni
Eigen::VectorXd System2D::calcChernNumbersFromBC(const Eigen::Vector2<std::size_t> &n, const H2D &dvxH, const H2D &dvyH)
{

    if (n.x() % 2 != 0 || n.y() % 2 != 0)
    {
        std::cerr << "Number of k-points must be even in each direction. Returning empty vector." << std::endl;
        return Eigen::VectorXd();
    }

    Eigen::VectorXd kx = Eigen::VectorXd::LinSpaced(n.x(), -M_PI, M_PI);
    Eigen::VectorXd ky = Eigen::VectorXd::LinSpaced(n.y(), -M_PI, M_PI);

    kx = kx.head(n.x() - 1).eval();
    ky = ky.head(n.y() - 1).eval();

    Eigen::VectorXd kx_t = kx.head(n.x() / 2 + 1);
    Eigen::VectorXd kx_b = kx.tail(n.x() / 2);

    Eigen::VectorXd ky_l = ky.head(n.y() / 2 + 1);
    Eigen::VectorXd ky_r = ky.tail(n.y() / 2);

    auto processGridSection = [this, dvxH, dvyH](const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec, Eigen::VectorXd &CN) -> void
    {
        CN.setZero();

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> SAES;

        Eigen::VectorXd BC(_n_bands); // Berry curvature

        Eigen::MatrixXcd dvxH_mat(_n_bands, _n_bands);
        Eigen::MatrixXcd dvyH_mat(_n_bands, _n_bands);

        Eigen::dcomplex v1x, v1y, v2x, v2y;

        for (auto kx : kx_vec)
        {
            for (auto ky : ky_vec)
            {
                SAES.compute(_H({kx, ky}, _p));

                dvxH_mat = dvxH({kx, ky}, _p);
                dvyH_mat = dvyH({kx, ky}, _p);

                BC.setZero();

                for (auto bi = 0; bi < _n_bands; ++bi) // bi - band index
                {

                    for (auto bj = 0; bj < _n_bands; ++bj)
                    {
                        if (bj == bi)
                            continue;

                        v1x = SAES.eigenvectors().col(bi).dot(dvxH_mat * SAES.eigenvectors().col(bj));
                        v1y = SAES.eigenvectors().col(bi).dot(dvyH_mat * SAES.eigenvectors().col(bj));

                        v2x = SAES.eigenvectors().col(bj).dot(dvxH_mat * SAES.eigenvectors().col(bi));
                        v2y = SAES.eigenvectors().col(bj).dot(dvyH_mat * SAES.eigenvectors().col(bi));

                        BC(bi) -= std::imag(v1x * v2y - v2x * v1y) / (SAES.eigenvalues()(bi) - SAES.eigenvalues()(bj)) / (SAES.eigenvalues()(bi) - SAES.eigenvalues()(bj));
                    }
                }

                CN += BC;
            }
        }

        CN *= (kx_vec(1) - kx_vec(0)) * (ky_vec(1) - ky_vec(0)) / (2.0 * M_PI);
    };

    Eigen::VectorXd CN_tl(_n_bands);
    Eigen::VectorXd CN_tr(_n_bands);
    Eigen::VectorXd CN_bl(_n_bands);
    Eigen::VectorXd CN_br(_n_bands);

    std::thread tl(processGridSection, kx_t, ky_l, std::ref(CN_tl));
    std::thread tr(processGridSection, kx_t, ky_r, std::ref(CN_tr));
    std::thread bl(processGridSection, kx_b, ky_l, std::ref(CN_bl));
    // std::thread br(processGridSection, kx_b, ky_r, std::ref(CN_br));

    // bottom right using main thread
    processGridSection(kx_b, ky_r, CN_br);

    tl.join();
    tr.join();
    bl.join();
    // br.join();

    return CN_tl + CN_tr + CN_bl + CN_br;
}

double System2D::calcCNUsingWilsonLoop(std::size_t n_dense, std::size_t n_sparse, double k_val)
{
    Eigen::VectorXd kx(n_dense + 2 * n_sparse + 1);
    kx.segment(0, n_sparse + 1) = Eigen::VectorXd::LinSpaced(n_sparse + 1, -M_PI, -k_val);
    kx.segment(n_sparse, n_dense + 1) = Eigen::VectorXd::LinSpaced(n_dense + 1, -k_val, k_val);
    kx.segment(n_sparse + n_dense, n_sparse + 1) = Eigen::VectorXd::LinSpaced(n_sparse + 1, k_val, M_PI);
    kx(kx.size() - 1) = -M_PI;

    Eigen::VectorXd ky = kx;

    std::size_t rows = kx.size();
    std::size_t cols = ky.size();

    Eigen::VectorX<Eigen::MatrixXcd> i_evecs(cols);  // Current row eigenvectors
    Eigen::VectorX<Eigen::MatrixXcd> ni_evecs(cols); // Next row eigenvectors

    Eigen::VectorX<Eigen::MatrixXcd> *i_evecs_ptr = &i_evecs;
    Eigen::VectorX<Eigen::MatrixXcd> *ni_evecs_ptr = &ni_evecs;

    // Compute eigenvectors for the first row
    for (auto j = 0; j < cols; ++j)
    {
        i_evecs(j) = _SAES.compute(_H({kx(0), ky(j)}, _p)).eigenvectors();
    }

    double CN = 0.0;
    std::size_t half_bands = _n_bands / 2;
    Eigen::MatrixXcd W(half_bands, half_bands);

    for (auto i = 0; i < rows - 1; ++i)
    {
        (*ni_evecs_ptr)(0) = _SAES.compute(_H({kx(i + 1), ky(0)}, _p)).eigenvectors();

        for (auto j = 0; j < cols - 1; ++j)
        {
            (*ni_evecs_ptr)(j + 1) = _SAES.compute(_H({kx(i + 1), ky(j + 1)}, _p)).eigenvectors(); // contains eigenvectors and eigenvalues for next i, next j

            W = (*i_evecs_ptr)(j).leftCols(half_bands).adjoint() * (*i_evecs_ptr)(j + 1).leftCols(half_bands);
            W *= (*i_evecs_ptr)(j + 1).leftCols(half_bands).adjoint() * (*ni_evecs_ptr)(j + 1).leftCols(half_bands);
            W *= (*ni_evecs_ptr)(j + 1).leftCols(half_bands).adjoint() * (*ni_evecs_ptr)(j).leftCols(half_bands);
            W *= (*ni_evecs_ptr)(j).leftCols(half_bands).adjoint() * (*i_evecs_ptr)(j).leftCols(half_bands);

            CN += std::arg(W.determinant());
        }

        std::swap(i_evecs_ptr, ni_evecs_ptr);
    }

    return CN / (2.0 * M_PI);
}

// tak sie nie da, bo nie ma sensu
double System2D::calcCNDiscreteky(std::size_t n_dense, std::size_t n_sparse, double k_val, std::size_t n_ky)
{
    Eigen::VectorXd kx(n_dense + 2 * n_sparse + 1);
    kx.segment(0, n_sparse + 1) = Eigen::VectorXd::LinSpaced(n_sparse + 1, -M_PI, -k_val);
    kx.segment(n_sparse, n_dense + 1) = Eigen::VectorXd::LinSpaced(n_dense + 1, -k_val, k_val);
    kx.segment(n_sparse + n_dense, n_sparse + 1) = Eigen::VectorXd::LinSpaced(n_sparse + 1, k_val, M_PI);
    kx(kx.size() - 1) = -M_PI;

    Eigen::MatrixXcd Ui = _SAES.compute(H_discrete_ky(kx(0), n_ky)).eigenvectors();
    Eigen::MatrixXcd Uip;

    Eigen::MatrixXcd *Ui_ptr = &Ui;
    Eigen::MatrixXcd *Uip_ptr = &Uip;

    std::size_t half_bands = _n_bands * n_ky / 2;

    Eigen::MatrixXcd W = Eigen::MatrixXcd::Identity(half_bands, half_bands);

    for (auto i = 0; i < kx.size() - 1; ++i)
    {
        *Uip_ptr = _SAES.compute(H_discrete_ky(kx(i + 1), n_ky)).eigenvectors();

        W *= (*Ui_ptr).leftCols(half_bands).adjoint() * (*Uip_ptr).leftCols(half_bands);

        std::swap(Ui_ptr, Uip_ptr);
    }

    double CN = std::arg(W.determinant()) / (2.0 * M_PI);

    return CN < -1e-3 ? CN + 1.0 : CN;
}

Eigen::MatrixXcd System2D::H_discrete_ky(double kx, std::size_t n_ky)
{
    Eigen::MatrixXcd H_ky = Eigen::MatrixXcd::Zero(_n_bands * n_ky, _n_bands * n_ky);

    Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_ky + 1, -M_PI, M_PI);
    ky_vec = ky_vec.head(n_ky).eval();

    for (int i = 0; i < n_ky; ++i)
        for (int j = i - 1; j < i + 2; ++j)
        {
            if (j < 0 || j >= n_ky)
                continue;

            for (double ky : ky_vec)
                H_ky.block(i * _n_bands, j * _n_bands, _n_bands, _n_bands) += std::exp(Eigen::dcomplex(0.0, ky * (i - j))) * _H({kx, ky}, _p);
        }

    H_ky /= n_ky;

    return H_ky;
}

// arma::sp_cx_mat System2D::armaH_discrete_ky(double kx, std::size_t n_ky)
// {
// arma::sp_cx_mat H(_n_bands * n_ky, _n_bands * n_ky);

// Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_ky + 1, -M_PI, M_PI);
// ky_vec = ky_vec.head(n_ky).eval();

// Eigen::MatrixXcd H_mat(_n_bands, _n_bands);
// arma::cx_mat armaH_mat(H_mat.data(), H_mat.rows(), H_mat.cols(), false);

// for (int i = 0; i < n_ky; ++i)
//     for (int j = i - 1; j <= i + 1; ++j)
//     {
//         if (j < 0 || j >= n_ky)
//             continue;

//         H_mat.setZero();

//         for (double ky : ky_vec)
//             H_mat += std::exp(Eigen::dcomplex(0.0, ky * (i - j))) * _H({kx, ky}, _p);

//         H.submat(i * _n_bands, j * _n_bands, arma::SizeMat(_n_bands, _n_bands)) = armaH_mat;
//     }

// return H / n_ky;
// }

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> System2D::H_discrete_ky_eigendecomposiiton(double kx, std::size_t n_ky, bool eigenvalues_olny, std::size_t n_eigs, double sigma)
{
    Eigen::MatrixXcd EigenH = H_discrete_ky(kx, n_ky);
    arma::cx_mat armaH(EigenH.data(), EigenH.rows(), EigenH.cols(), false);
    arma::sp_cx_mat armaHsp(armaH);

    arma::eigs_opts opts;
    opts.tol = meV2au(1e-6);
    opts.subdim = 3 * n_eigs + 10;

    arma::cx_vec eigval(n_eigs);
    Eigen::VectorXd evals(n_eigs);

    if (eigenvalues_olny)
    {
        arma::eigs_gen(eigval, armaHsp, n_eigs, sigma, opts);

        for (auto i = 0; i < n_eigs; ++i)
        {
            evals(i) = eigval(i).real();
        }

        return {evals, Eigen::MatrixXcd()};
    }
    else
    {
        arma::cx_mat eigvec(_n_bands * n_ky, n_eigs);
        arma::eigs_gen(eigval, eigvec, armaHsp, n_eigs, sigma, opts);

        Eigen::MatrixXcd evecs(n_ky * _n_bands, n_eigs);

        for (auto i = 0; i < n_eigs; ++i)
        {
            evals(i) = eigval(i).real();

            for (auto j = 0; j < n_ky * _n_bands; ++j)
            {
                evecs(j, i) = eigvec(j, i);
            }
        }

        return {evals, evecs};
    }
}

Eigen::MatrixXcd System2D::H_discrete(std::size_t n_kx, std::size_t n_ky)
{
    Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(_n_bands * n_kx * n_ky, _n_bands * n_kx * n_ky);

    Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_kx + 1, -M_PI, M_PI);
    Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_ky + 1, -M_PI, M_PI);

    kx_vec = kx_vec.head(n_kx).eval();
    ky_vec = ky_vec.head(n_ky).eval();

    for (int i = 0; i < n_kx; ++i)
        for (int j = i - 1; j <= i + 1; ++j)
        {
            if (j < 0 || j >= n_kx)
                continue;

            std::size_t ix = i * n_ky * _n_bands;
            std::size_t jx = j * n_ky * _n_bands;

            for (int k = 0; k < n_ky; ++k)
                for (int l = k - 1; l <= k + 1; ++l)
                {
                    if (l < 0 || l >= n_ky)
                        continue;

                    for (double kx : kx_vec)
                        for (double ky : ky_vec)
                        {
                            double phase = kx * (i - j) + ky * (k - l);
                            H.block(ix + k * _n_bands, jx + l * _n_bands, _n_bands, _n_bands) += std::exp(1i * phase) * _H({kx, ky}, _p);
                        }
                }
        }

    return H / (n_kx * n_ky);
}

// arma::sp_cx_mat System2D::armaH_discrete(std::size_t n_kx, std::size_t n_ky)
// {
//     arma::sp_cx_mat H(_n_bands * n_kx * n_ky, _n_bands * n_kx * n_ky);

//     Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_kx + 1, -M_PI, M_PI);
//     Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_ky + 1, -M_PI, M_PI);

//     kx_vec = kx_vec.head(n_kx).eval();
//     ky_vec = ky_vec.head(n_ky).eval();

//     Eigen::MatrixXcd H_mat(_n_bands, _n_bands);
//     arma::cx_mat armaH_mat(H_mat.data(), H_mat.rows(), H_mat.cols(), false);

//     for (int i = 0; i < n_kx; ++i)
//         for (int j = i - 1; j <= i + 1; ++j)
//         {
//             if (j < 0 || j >= n_kx)
//                 continue;

//             std::size_t ix = i * n_ky * _n_bands;
//             std::size_t jx = j * n_ky * _n_bands;

//             for (int k = 0; k < n_ky; ++k)
//                 for (int l = k - 1; l <= k + 1; ++l)
//                 {
//                     if (l < 0 || l >= n_ky)
//                         continue;

//                     H_mat.setZero();

//                     for (double kx : kx_vec)
//                         for (double ky : ky_vec)
//                             H_mat += std::exp(Eigen::dcomplex(0.0, kx * (i - j) + ky * (k - l))) * _H({kx, ky}, _p);

//                     H.submat(ix + k * _n_bands, jx + l * _n_bands, arma::SizeMat(_n_bands, _n_bands)) = armaH_mat;
//                 }
//         }

//     return H / (n_kx * n_ky);
// }

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> System2D::H_discrete_eigendecomposiiton(std::size_t n_kx, std::size_t n_ky, bool eigenvalues_olny, std::size_t n_eigs, double sigma)
{
    Eigen::MatrixXcd EigenH = H_discrete(n_kx, n_ky);
    arma::cx_mat armaH(EigenH.data(), EigenH.rows(), EigenH.cols(), false);
    arma::sp_cx_mat armaHsp(armaH);

    arma::eigs_opts opts;
    opts.tol = meV2au(1e-6);
    opts.subdim = 3 * n_eigs + 10;

    arma::cx_vec eigval;
    Eigen::VectorXd evals(n_eigs);

    if (eigenvalues_olny)
    {
        arma::eigs_gen(eigval, armaHsp, n_eigs, sigma, opts);

        for (auto i = 0; i < n_eigs; ++i)
        {
            evals(i) = eigval(i).real();
        }

        return {evals, Eigen::MatrixXcd()};
    }
    else
    {
        arma::cx_mat eigvec;
        arma::eigs_gen(eigval, eigvec, armaHsp, n_eigs, sigma, opts);

        Eigen::MatrixXcd evecs(n_kx * n_ky * _n_bands, n_eigs);

        for (auto i = 0; i < n_eigs; ++i)
        {
            evals(i) = eigval(i).real();

            for (auto j = 0; j < n_kx * n_ky * _n_bands; ++j)
            {
                evecs(j, i) = eigvec(j, i);
            }
        }

        return {evals, evecs};
    }
}
