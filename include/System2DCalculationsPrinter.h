#ifndef SYSTEM2DCALCULATIONSPRINTER_H
#define SYSTEM2DCALCULATIONSPRINTER_H

#include "System2DCalculations.h"

class System2DCalculationsPrinter
{
public:
    System2DCalculationsPrinter() = delete;
    System2DCalculationsPrinter(const System2DCalculations &calc) : _calc(calc) {}

    void printBandStructure(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec);
    void printBandStructureSlice(const std::string &output_filename, const Eigen::VectorXd &k_vec, int axis, double k0 = 0.0);
    void printBandStructureSlice_normal(const std::string &output_filename, const Eigen::VectorXd &k_vec, int axis, double k0 = 0.0);
    void printBandStructure_discrete_kx(const std::string &output_filename, const Eigen::VectorXd &ky_vec, std::size_t n_kx);
    void printBandStructure_discrete_kx_normal(const std::string &output_filename, const Eigen::VectorXd &ky_vec, std::size_t n_kx);
    void printBandStructure_discrete_ky(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky);
    void printBandStructure_discrete_ky_normal(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky);
    void printBandStructure_sparse_discrete_kx(const std::string &output_filename, const Eigen::VectorXd &ky_vec, std::size_t n_kx, std::size_t n_eigs = 30, double sigma = 0.0);
    void printBandStructure_sparse_discrete_kx_normal(const std::string &output_filename, const Eigen::VectorXd &ky_vec, std::size_t n_kx, std::size_t n_eigs = 30, double sigma = 0.0);
    void printBandStructure_sparse_discrete_ky(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky, std::size_t n_eigs = 30, double sigma = 0.0);
    void printBandStructure_sparse_discrete_ky_normal(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky, std::size_t n_eigs = 30, double sigma = 0.0);

    void printBandStructure_orbital_type(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec);
    void printBandStructureSlice_orbital_type(const std::string &output_filename, const Eigen::VectorXd &k_vec, int axis, double k0 = 0.0);
    void printBandStructureSlice_normal_orbital_type(const std::string &output_filename, const Eigen::VectorXd &k_vec, int axis, double k0 = 0.0);
    void printBandStructure_discrete_kx_orbital_type(const std::string &output_filename, const Eigen::VectorXd &ky_vec, std::size_t n_kx);
    void printBandStructure_discrete_kx_normal_orbital_type(const std::string &output_filename, const Eigen::VectorXd &ky_vec, std::size_t n_kx);
    void printBandStructure_discrete_ky_orbital_type(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky);
    void printBandStructure_discrete_ky_normal_orbital_type(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky);
    void printBandStructure_sparse_discrete_kx_orbital_type(const std::string &output_filename, const Eigen::VectorXd &ky_vec, std::size_t n_kx, std::size_t n_eigs = 30, double sigma = 0.0);
    void printBandStructure_sparse_discrete_kx_normal_orbital_type(const std::string &output_filename, const Eigen::VectorXd &ky_vec, std::size_t n_kx, std::size_t n_eigs = 30, double sigma = 0.0);
    void printBandStructure_sparse_discrete_ky_orbital_type(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky, std::size_t n_eigs = 30, double sigma = 0.0);
    void printBandStructure_sparse_discrete_ky_normal_orbital_type(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky, std::size_t n_eigs = 30, double sigma = 0.0);

    void printProbDen_sparse_discrete(const std::string &output_filename, std::size_t nk_x, std::size_t nk_y, double E);

    void printAbsDelta(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec);
    void printAbsDeltaAlongContour(const std::string &output_filename, const std::vector<Point2D> &contour);

    void printAbelianBerryCurvature(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec);
    void printBerryCurvatureFromWilsonLoop(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec);

    void printWilsonLoopSpectrum(const std::string &output_filename, std::size_t nk, std::size_t nk_loop, double k_max = M_PI);

private:
    System2DCalculations _calc;
};

#endif // SYSTEM2DCALCULATIONS_H