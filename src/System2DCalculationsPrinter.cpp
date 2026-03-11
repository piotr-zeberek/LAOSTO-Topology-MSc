#include "System2DCalculationsPrinter.h"
#include "utils.h"

#include <iostream>
#include <iomanip>
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

void System2DCalculationsPrinter::printBandStructureSlice_normal(const std::string &output_filename, const Eigen::VectorXd &k_vec, int axis, double k0)
{
    std::ofstream output_file(output_filename);

    if (axis == 0)
    {
        for (auto kx : k_vec)
            output_file << kx << " " << _calc.eigenvals_normal(kx, k0).transpose() / meV2au(1.0) << "\n";
    }
    else if (axis == 1)
    {
        for (auto ky : k_vec)
            output_file << ky << " " << _calc.eigenvals_normal(k0, ky).transpose() / meV2au(1.0) << "\n";
    }
    else
    {
        std::cerr << "Invalid axis. Choose 0 or 1." << std::endl;
    }
}

void System2DCalculationsPrinter::printBandStructure_discrete_kx(const std::string &output_filename, const Eigen::VectorXd &ky_vec, std::size_t n_kx)
{
    std::ofstream output_file(output_filename);

    for (auto ky : ky_vec)
    {
        auto evals = _calc.eigenvals_discrete_kx(n_kx, ky);

        output_file << ky << " " << evals.transpose() / meV2au(1.0) << std::endl;
    }
}

void System2DCalculationsPrinter::printBandStructure_discrete_kx_normal(const std::string &output_filename, const Eigen::VectorXd &ky_vec, std::size_t n_kx)
{
    std::ofstream output_file(output_filename);

    for (auto ky : ky_vec)
    {
        auto evals = _calc.eigenvals_discrete_kx_normal(n_kx, ky);

        output_file << ky << " " << evals.transpose() / meV2au(1.0) << std::endl;
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

void System2DCalculationsPrinter::printBandStructure_discrete_ky_normal(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
    {
        auto evals = _calc.eigenvals_discrete_ky_normal(kx, n_ky);

        output_file << kx << " " << evals.transpose() / meV2au(1.0) << std::endl;
    }
}

void System2DCalculationsPrinter::printBandStructure_sparse_discrete_kx(const std::string &output_filename, const Eigen::VectorXd &ky_vec, std::size_t n_kx, std::size_t n_eigs, double sigma)
{
    std::ofstream output_file(output_filename);

    for (auto ky : ky_vec)
    {
        auto evals = _calc.eigenvals_sparse_discrete_kx(n_kx, ky, n_eigs, sigma);

        output_file << ky << " " << evals.transpose() / meV2au(1.0) << std::endl;
    }
}

void System2DCalculationsPrinter::printBandStructure_sparse_discrete_kx_normal(const std::string &output_filename, const Eigen::VectorXd &ky_vec, std::size_t n_kx, std::size_t n_eigs, double sigma)
{
    std::ofstream output_file(output_filename);

    for (auto ky : ky_vec)
    {
        auto evals = _calc.eigenvals_sparse_discrete_kx_normal(n_kx, ky, n_eigs, sigma);

        output_file << ky << " " << evals.transpose() / meV2au(1.0) << std::endl;
    }
}

void System2DCalculationsPrinter::printBandStructure_sparse_discrete_ky(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky, std::size_t n_eigs, double sigma)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
    {
        auto evals = _calc.eigenvals_sparse_discrete_ky(kx, n_ky, n_eigs, sigma);

        output_file << kx << " " << evals.transpose() / meV2au(1.0) << std::endl;
    }
}

void System2DCalculationsPrinter::printBandStructure_sparse_discrete_ky_normal(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky, std::size_t n_eigs, double sigma)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
    {
        auto evals = _calc.eigenvals_sparse_discrete_ky_normal(kx, n_ky, n_eigs, sigma);

        output_file << kx << " " << evals.transpose() / meV2au(1.0) << std::endl;
    }
}

void System2DCalculationsPrinter::printBandStructure_orbital_type(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
        for (auto ky : ky_vec)
        {
            auto [evals, evecs] = _calc.eigen(kx, ky);

            output_file << kx << " " << ky;

            for (auto i = 0; i < evals.size(); ++i)
            {
                output_file << " " << evals(i) / meV2au(1.0) << " " << orbital_prob_den(evecs.col(i), _calc.system().n_bands * 2).transpose();
            }
            output_file << "\n";
        }
}
void System2DCalculationsPrinter::printBandStructureSlice_orbital_type(const std::string &output_filename, const Eigen::VectorXd &k_vec, int axis, double k0)
{
    std::ofstream output_file(output_filename);

    for (auto k : k_vec)
    {
        double kx = (axis == 0) ? k : k0;
        double ky = (axis == 1) ? k : k0;

        auto [evals, evecs] = _calc.eigen(kx, ky);
        output_file << k;
        for (auto i = 0; i < evals.size(); ++i)
        {
            output_file << std::setprecision(9) << " " << evals(i) / meV2au(1.0) << " " << orbital_prob_den(evecs.col(i), _calc.system().n_bands * 2).transpose();
        }
        output_file << "\n";
    }
}
void System2DCalculationsPrinter::printBandStructureSlice_normal_orbital_type(const std::string &output_filename, const Eigen::VectorXd &k_vec, int axis, double k0)
{
    std::ofstream output_file(output_filename);

    for (auto k : k_vec)
    {
        double kx = (axis == 0) ? k : k0;
        double ky = (axis == 1) ? k : k0;

        auto [evals, evecs] = _calc.eigen_normal(kx, ky);
        output_file << k;
        for (auto i = 0; i < evals.size(); ++i)
        {
            output_file << " " << evals(i) / meV2au(1.0) << " " << orbital_prob_den(evecs.col(i), _calc.system().n_bands).transpose();
        }
        output_file << "\n";
    }
}

void System2DCalculationsPrinter::printBandStructure_discrete_kx_orbital_type(const std::string &output_filename, const Eigen::VectorXd &ky_vec, std::size_t n_kx)
{
    std::ofstream output_file(output_filename);

    for (auto ky : ky_vec)
    {
        auto [evals, evecs] = _calc.eigen_discrete_kx(n_kx, ky);
        output_file << ky;
        for (auto i = 0; i < evals.size(); ++i)
        {
            output_file << " " << evals(i) / meV2au(1.0) << " " << orbital_prob_den(evecs.col(i), _calc.system().n_bands * 2).transpose();
        }
        output_file << "\n";
    }
}

void System2DCalculationsPrinter::printBandStructure_discrete_kx_normal_orbital_type(const std::string &output_filename, const Eigen::VectorXd &ky_vec, std::size_t n_kx)
{
    std::ofstream output_file(output_filename);

    for (auto ky : ky_vec)
    {
        auto [evals, evecs] = _calc.eigen_discrete_kx_normal(n_kx, ky);
        output_file << ky;
        for (auto i = 0; i < evals.size(); ++i)
        {
            output_file << " " << evals(i) / meV2au(1.0) << " " << orbital_prob_den(evecs.col(i), _calc.system().n_bands).transpose();
        }
        output_file << "\n";
    }
}

void System2DCalculationsPrinter::printBandStructure_discrete_ky_orbital_type(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
    {
        auto [evals, evecs] = _calc.eigen_discrete_ky(kx, n_ky);
        output_file << kx;
        for (auto i = 0; i < evals.size(); ++i)
        {
            output_file << " " << evals(i) / meV2au(1.0) << " " << orbital_prob_den(evecs.col(i), _calc.system().n_bands * 2).transpose();
        }
        output_file << "\n";
    }
}
void System2DCalculationsPrinter::printBandStructure_discrete_ky_normal_orbital_type(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
    {
        auto [evals, evecs] = _calc.eigen_discrete_ky_normal(kx, n_ky);
        output_file << kx;
        for (auto i = 0; i < evals.size(); ++i)
        {

            output_file << " " << evals(i) / meV2au(1.0) << " " << orbital_prob_den(evecs.col(i), _calc.system().n_bands).transpose();
        }
        output_file << "\n";
    }
}

void System2DCalculationsPrinter::printBandStructure_sparse_discrete_kx_orbital_type(const std::string &output_filename, const Eigen::VectorXd &ky_vec, std::size_t n_kx, std::size_t n_eigs, double sigma)
{
    std::ofstream output_file(output_filename);

    for (auto ky : ky_vec)
    {
        auto [evals, evecs] = _calc.eigen_sparse_discrete_kx(n_kx, ky, n_eigs, sigma);
        output_file << ky;
        for (auto i = 0; i < evals.size(); ++i)
        {
            output_file << " " << evals(i) / meV2au(1.0) << " " << orbital_prob_den(evecs.col(i), _calc.system().n_bands * 2).transpose();
        }
        output_file << "\n";
    }
}
void System2DCalculationsPrinter::printBandStructure_sparse_discrete_kx_normal_orbital_type(const std::string &output_filename, const Eigen::VectorXd &ky_vec, std::size_t n_kx, std::size_t n_eigs, double sigma)
{
    std::ofstream output_file(output_filename);

    for (auto ky : ky_vec)
    {
        auto [evals, evecs] = _calc.eigen_sparse_discrete_kx_normal(n_kx, ky, n_eigs, sigma);
        output_file << ky;
        for (auto i = 0; i < evals.size(); ++i)
        {
            output_file << " " << evals(i) / meV2au(1.0) << " " << orbital_prob_den(evecs.col(i), _calc.system().n_bands).transpose();
        }
        output_file << "\n";
    }
}
void System2DCalculationsPrinter::printBandStructure_sparse_discrete_ky_orbital_type(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky, std::size_t n_eigs, double sigma)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
    {
        auto [evals, evecs] = _calc.eigen_sparse_discrete_ky(kx, n_ky, n_eigs, sigma);
        output_file << kx;
        for (auto i = 0; i < evals.size(); ++i)
        {
            output_file << " " << evals(i) / meV2au(1.0) << " " << orbital_prob_den(evecs.col(i), _calc.system().n_bands * 2).transpose();
        }
        output_file << std::endl;
    }
}

void System2DCalculationsPrinter::printBandStructure_sparse_discrete_ky_normal_orbital_type(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky, std::size_t n_eigs, double sigma)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
    {
        auto [evals, evecs] = _calc.eigen_sparse_discrete_ky_normal(kx, n_ky, n_eigs, sigma);
        output_file << kx;
        for (auto i = 0; i < evals.size(); ++i)
        {
            output_file << " " << evals(i) / meV2au(1.0) << " " << orbital_prob_den(evecs.col(i), _calc.system().n_bands).transpose();
        }
        output_file << "\n";
    }
}

void System2DCalculationsPrinter::printProbDen_sparse_discrete(const std::string &output_filename, std::size_t nk_x, std::size_t nk_y, double E)
{
    std::ofstream output_file(output_filename);

    auto [evals, evecs] = _calc.eigen_sparse_discrete(nk_x, nk_y, 1, E);
    Eigen::VectorXd prob_den_orbitals = evecs.col(0).cwiseAbs2();

    output_file << "# E = " << E / meV2au(1) << " eval = " << evals(0) / meV2au(1) << std::endl;

    for (auto i = 0; i < nk_x; ++i)
    {
        for (auto j = 0; j < nk_y; ++j)
        {
            output_file << i << " " << j << " " << prob_den_orbitals.segment((i * nk_y + j) * 2 * _calc.system().n_bands, 2 * _calc.system().n_bands).sum() << std::endl;
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

void System2DCalculationsPrinter::printWilsonLoopSpectrum(const std::string &output_filename, std::size_t nk, std::size_t nk_loop, double k_max)
{
    Eigen::VectorXd k_vec = Eigen::VectorXd::LinSpaced(nk, -k_max, k_max);

    std::ofstream output_file(output_filename);

    for (auto k : k_vec)
    {
        output_file << k << " " << _calc.WilsonLoopSpectrum(nk_loop, 0, k).transpose() << " " << _calc.WilsonLoopSpectrum(nk_loop, 1, k).transpose() << std::endl;
    }
}