#include "System2DCalculationsPrinter.h"
#include "utils.h"

#include <iostream>
#include <fstream>

void System2DCalculationsPrinter::printBandStructure(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
        for (auto ky : ky_vec)
            output_file << kx << " " << ky << " " << _calc.eigenvals(kx, ky).transpose() / meV2au(1.0) << "\n";
}

void System2DCalculationsPrinter::printBandStructureSlice(const std::string &output_filename, const Eigen::VectorXd &k_vec, int axis, double k0)
{
    std::ofstream output_file(output_filename);

    if (axis == 0)
    {
        for (auto kx : k_vec)
            output_file << kx << " " << _calc.eigenvals(kx, k0).transpose() / meV2au(1.0) << "\n";
    }
    else if (axis == 1)
    {
        for (auto ky : k_vec)
            output_file << ky << " " << _calc.eigenvals(k0, ky).transpose() / meV2au(1.0) << "\n";
    }
    else
    {
        std::cerr << "Invalid axis. Choose 0 or 1." << std::endl;
    }
}

void System2DCalculationsPrinter::printBandStructure_discrete_ky(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
    {
        auto evals = _calc.eigenvals_discrete_ky(kx, n_ky);

        output_file << kx << " " << evals.transpose() / meV2au(1.0) << std::endl;
    }
}

void System2DCalculationsPrinter::printBandStructure_sparse_discrete_ky(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky, std::size_t n_eigs, double sigma)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
    {
        auto evals = _calc.eigenvals_sprase_discrete_ky(kx, n_ky, n_eigs, sigma);

        output_file << kx << " " << evals.transpose() / meV2au(1.0) << std::endl;
    }
}

void System2DCalculationsPrinter::printProbDen_sparse_discrete(const std::string &output_filename, std::size_t nk_x, std::size_t nk_y, double E)
{
    std::ofstream output_file(output_filename);

    auto [evals, evecs] = _calc.eigen_sprase_discrete(nk_x, nk_y, 1, E);
    Eigen::VectorXd prob_den_orbitals = evecs.col(0).cwiseAbs2();

    output_file << "# E = " << E / meV2au(1) << " eval = " << evals(0) / meV2au(1) << std::endl;

    for (auto i = 0; i < nk_x; ++i)
    {
        for (auto j = 0; j < nk_y; ++j)
        {
            output_file << i << " " << j << " " << prob_den_orbitals.segment((i * nk_y + j) * 2*_calc.system().n_bands, 2*_calc.system().n_bands).sum() << std::endl;
        }
    }
}

void System2DCalculationsPrinter::printAbsDelta(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
        for (auto ky : ky_vec)
            output_file << kx << " " << ky << " " << _calc.AbsDelta(kx, ky).transpose() / meV2au(1) << std::endl;
}

void System2DCalculationsPrinter::printAbsDeltaAlongContour(const std::string &output_filename, const std::vector<Point2D> &contour)
{
    std::ofstream output_file(output_filename);

    for (auto k : contour)
    {
        output_file << k.x() << " " << k.y() << " " << _calc.AbsDelta(k.x(), k.y()).transpose() / meV2au(1) << "\n";
    }
}

// // bujda na resorach te wyniki - fazy tancza w kazda strone
// void System2DCalculationsPrinter::printDeltaFromUnitaryTransformation(const std::string &delta_filename, const std::string &DT_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec)
// {
//     std::ofstream delta_file(delta_filename);
//     std::ofstream DT_file(DT_filename);

//     std::size_t half_bands = _n_bands / 2;
//     std::size_t orbital_bands = half_bands / 2;

//     Eigen::MatrixXcd H(_n_bands, _n_bands);
//     Eigen::MatrixXcd U = Eigen::MatrixXcd::Zero(_n_bands, _n_bands);

//     Eigen::MatrixXcd delta_matrix(half_bands, half_bands);
//     Eigen::VectorXcd delta_matrix_vector_view(half_bands * half_bands);

//     for (auto kx : kx_vec)
//     {
//         for (auto ky : ky_vec)
//         {
//             H = _sys.H({kx, ky}, _p);
//             U.topLeftCorner(half_bands, half_bands) = _SAES.compute(H.topLeftCorner(half_bands, half_bands)).eigenvectors();
//             U.bottomRightCorner(half_bands, half_bands) = _SAES.compute(H.bottomRightCorner(half_bands, half_bands)).eigenvectors().rowwise().reverse();

//             // //stare proby cechowania fazy
//             // std::size_t i = 0;
//             // std::size_t j = 0;
//             // U.cwiseAbs2().maxCoeff(&i, &j);
//             // U /= U(i, j);
//             // U.colwise().normalize();

//             delta_matrix = (U.adjoint() * H * U).topRightCorner(half_bands, half_bands);

//             delta_file << kx << " " << ky << " ";

//             for (auto i = 0; i < half_bands; ++i)
//             {
//                 for (auto j = 0; j < half_bands; ++j)
//                 {
//                     delta_file << std::abs(delta_matrix(i, j)) / meV2au(1) << " " << std::arg(delta_matrix(i, j)) << " ";
//                 }
//             }

//             delta_file << std::endl;

//             Eigen::dcomplex det = delta_matrix.determinant();
//             Eigen::dcomplex tr = delta_matrix.trace();

//             DT_file << kx << " " << ky << " "
//                     << std::abs(det) / std::pow(meV2au(1), half_bands) << " " << std::arg(det) << " "
//                     << std::abs(tr) / meV2au(1) << " " << std::arg(tr) << std::endl;
//         }
//     }
// }

void System2DCalculationsPrinter::printAbelianBerryCurvature(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
    {
        for (auto ky : ky_vec)
        {
            output_file << kx << " " << ky << " " << _calc.AbelianBerryCurvature(kx, ky).transpose() << "\n";
        }
    }
}

void System2DCalculationsPrinter::printBerryCurvatureFromWilsonLoop(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
    {
        for (auto ky : ky_vec)
        {
            output_file << kx << " " << ky << " " << _calc.BerryCurvatureFromWilsonLoop(kx, ky) << "\n";
        }
    }
}

// void System2DCalculationsPrinter::printMatrixBerryCurvatureTrace(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec, double BC_trace_skip)
// {
//     std::ofstream output_file(output_filename);

//     double BC_trace;

//     for (auto i = 0; i < kx_vec.size(); ++i)
//     {
//         for (auto j = 0; j < ky_vec.size(); ++j)
//         {
//             BC_trace = calcMatrixBerryCurvatureTrace({kx_vec(i), ky_vec(j)});

//             // if (std::abs(BC_trace) > BC_trace_skip)
//             // {
//             //     BC_trace = 0.0;
//             // }

//             output_file << kx_vec(i) << " " << ky_vec(j) << " " << BC_trace << "\n";
//         }
//     }
// }

void System2DCalculationsPrinter::printWilsonLoopSpectrum(const std::string &output_filename, std::size_t nk, std::size_t nk_loop, double k_max)
{
    Eigen::VectorXd k_vec = Eigen::VectorXd::LinSpaced(nk, -k_max, k_max);

    std::ofstream output_file(output_filename);

    for (auto k : k_vec)
    {
        output_file << k << " " << _calc.WilsonLoopSpectrum(nk_loop, 0, k).transpose() << " " << _calc.WilsonLoopSpectrum(nk_loop, 1, k).transpose() << std::endl;
    }
}