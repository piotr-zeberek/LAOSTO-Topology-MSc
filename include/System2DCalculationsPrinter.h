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
    void printBandStructure_discrete_ky(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky);
    void printBandStructure_sparse_discrete_ky(const std::string &output_filename, const Eigen::VectorXd &kx_vec, std::size_t n_ky, std::size_t n_eigs = 30, double sigma = 0.0);
    
    void printProbDen_sparse_discrete(const std::string &output_filename, std::size_t nk_x, std::size_t nk_y, double E);

    void printAbsDelta(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec);
    void printAbsDeltaAlongContour(const std::string &output_filename, const std::vector<Point2D> &contour);
    void printDeltaFromUnitaryTransformation(const std::string &delta_filename, const std::string &DT_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec);

    void printAbelianBerryCurvature(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec);
    void printBerryCurvatureFromWilsonLoop(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec);
    void printMatrixBerryCurvatureTrace(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec, double BC_trace_skip = 1e3);

    void printWilsonLoopSpectrum(const std::string &output_filename, std::size_t nk, std::size_t nk_loop, double k_max = M_PI);

private:
    System2DCalculations _calc;
};

#endif // SYSTEM2DCALCULATIONS_H