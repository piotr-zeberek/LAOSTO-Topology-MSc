#ifndef SYSTEM2DCALCULATIONS_H
#define SYSTEM2DCALCULATIONS_H

#include "System2D.h"

#include <Eigen/SVD>

using Point2D = Eigen::Vector2d;

class System2DCalculations
{
public:
    System2DCalculations() = delete;
    System2DCalculations(System2D &sys) : _sys(sys), _SAES(sys.n_bands * 2), _SVD(sys.n_bands * 2, sys.n_bands * 2) {}

    const System2D &system() { return _sys; }

    Hamiltonian dHdkx(double kx, double ky, double dk = 1e-6) const;
    Hamiltonian dHdky(double kx, double ky, double dk = 1e-6) const;
    Hamiltonian dHdkx_normal(double kx, double ky, double dk = 1e-6) const;
    Hamiltonian dHdky_normal(double kx, double ky, double dk = 1e-6) const;
    Eigen::VectorXd generate_k_vec(std::size_t n_dense, std::size_t n_sparse, double k_val) const;

    Eigen::VectorXd eigenvals(double kx, double ky);
    Eigen::VectorXd eigenvals_normal(double kx, double ky);
    Eigen::VectorXd eigenvals_discrete_ky(double kx, std::size_t n_ky);
    Eigen::VectorXd eigenvals_discrete_ky_normal(double kx, std::size_t n_ky);
    Eigen::VectorXd eigenvals_discrete(std::size_t n_kx, std::size_t n_ky);
    Eigen::VectorXd eigenvals_discrete_normal(std::size_t n_kx, std::size_t n_ky);
    
    Eigen::VectorXd eigenvals_sparse_discrete_ky(double kx, std::size_t n_ky, std::size_t n_eigs = 30, double sigma = 0.0);
    Eigen::VectorXd eigenvals_sparse_discrete_ky_normal(double kx, std::size_t n_ky, std::size_t n_eigs = 30, double sigma = 0.0);
    Eigen::VectorXd eigenvals_sparse_discrete(std::size_t n_kx, std::size_t n_ky, std::size_t n_eigs = 30, double sigma = 0.0);
    Eigen::VectorXd eigenvals_sparse_discrete_normal(std::size_t n_kx, std::size_t n_ky, std::size_t n_eigs = 30, double sigma = 0.0);

    Eigen::MatrixXcd eigenvecs(double kx, double ky);
    Eigen::MatrixXcd eigenvecs_normal(double kx, double ky);
    Eigen::MatrixXcd eigenvecs_discrete_ky(double kx, std::size_t n_ky);
    Eigen::MatrixXcd eigenvecs_discrete_ky_normal(double kx, std::size_t n_ky);
    Eigen::MatrixXcd eigenvecs_discrete(std::size_t n_kx, std::size_t n_ky);
    Eigen::MatrixXcd eigenvecs_discrete_normal(std::size_t n_kx, std::size_t n_ky);

    std::pair<Eigen::VectorXd, Eigen::MatrixXcd> eigen(double kx, double ky);
    std::pair<Eigen::VectorXd, Eigen::MatrixXcd> eigen_normal(double kx, double ky);
    std::pair<Eigen::VectorXd, Eigen::MatrixXcd> eigen_discrete_ky(double kx, std::size_t n_ky);
    std::pair<Eigen::VectorXd, Eigen::MatrixXcd> eigen_discrete_ky_normal(double kx, std::size_t n_ky);
    std::pair<Eigen::VectorXd, Eigen::MatrixXcd> eigen_discrete(std::size_t n_kx, std::size_t n_ky);
    std::pair<Eigen::VectorXd, Eigen::MatrixXcd> eigen_discrete_normal(std::size_t n_kx, std::size_t n_ky);
    
    std::pair<Eigen::VectorXd, Eigen::MatrixXcd> eigen_sparse_discrete_ky(double kx, std::size_t n_ky, std::size_t n_eigs = 30, double sigma = 0.0);
    std::pair<Eigen::VectorXd, Eigen::MatrixXcd> eigen_sparse_discrete_ky_normal(double kx, std::size_t n_ky, std::size_t n_eigs = 30, double sigma = 0.0);
    std::pair<Eigen::VectorXd, Eigen::MatrixXcd> eigen_sparse_discrete(std::size_t n_kx, std::size_t n_ky, std::size_t n_eigs = 30, double sigma = 0.0);
    std::pair<Eigen::VectorXd, Eigen::MatrixXcd> eigen_sparse_discrete_normal(std::size_t n_kx, std::size_t n_ky, std::size_t n_eigs = 30, double sigma = 0.0);

    Eigen::VectorXd AbsDelta(double kx, double ky);

    Eigen::VectorXd AbelianBerryCurvature(double kx, double ky);
    double BerryCurvatureFromWilsonLoop(double kx, double ky, double dk = 1e-2);
    double MatrixBerryCurvatureTrace(double kx, double ky, double dk = 1e-6);

    Eigen::VectorXd WilsonLoopSpectrum(std::size_t n, int axis, double k0 = 0.0);

    Eigen::VectorXd ChernNumbersUsingAbelianBerryCurvature(std::size_t n_dense, std::size_t n_sparse, double k_val);
    double ChernNumberUsingBerryCurvatureFromWilsonLoop(std::size_t n_dense, std::size_t n_sparse, double k_val);
    double ChernNumberUsingWilsonLoop(std::size_t n_dense, std::size_t n_sparse, double k_val);
    // double ChernNumberUsingMatrixBerryCurvatureTrace(std::size_t n_dense, std::size_t n_sparse, double k_val, double CN_skip = 1e3);

    double ChernNumberUsingWilsonLoop_discrete_ky(std::size_t n_dense, std::size_t n_sparse, double k_val, std::size_t n_ky);

    std::vector<std::vector<Point2D>> FSContours(double E = 0.0, double dk = 1e-4, double eps = 1e-9, double kx_min = -M_PI, double kx_max = 0.0, std::size_t n_kx = 1001);

    // Eigen::VectorXd calcChernNumbersDenserCenter(std::size_t n_dense, std::size_t n_sparse, double k_val);
    // Eigen::VectorXd calcChernNumbersWithCustomGrid(const Eigen::VectorXd &kx, const Eigen::VectorXd &ky);

private:
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> _SAES;
    Eigen::JacobiSVD<Eigen::MatrixXcd, Eigen::ComputeFullU | Eigen::ComputeFullV> _SVD;
    System2D &_sys;
};

#endif // SYSTEM2DCALCULATIONS_H