#ifndef SYSTEM2D_H
#define SYSTEM2D_H

#include <functional>
#include <string>
#include <tuple>

#include <Eigen/Core>
#include <Eigen/Eigen>

struct Parameters
{
    // Hopping parameters: l - light, h - heavy, d - coupling between the dxz /dyz
    double tl{};
    double th{};
    double td{};

    // Energy difference between xy and xz/yz orbitals
    double delta_E{};

    // Atomic spin-orbit coupling
    double delta_SO{};

    // Rashba spin-orbit coupling
    double delta_RSO{};

    // Superconducting energy gap
    double delta_SC{};

    // Chemical potential
    double mu{};

    // Lande g-factor
    double g_Lande{};

    // External magnetic field
    double Bx{};
    double By{};
    double Bz{};

    // effective mass
    double m{};

    // Rashba coupling
    double alpha{};

    // Hopping amplitude
    double t{};

    // lattice constant
    double a{};
};

using H2D = std::function<Eigen::MatrixXcd(const Eigen::Vector2d &k, const Parameters &p)>;
using Point2D = Eigen::Vector2d;

class System2D
{
public:
    System2D() = delete;
    System2D(const H2D &H, const Parameters &p = Parameters());

    void setHamiltonian(const H2D &H)
    {
        _H = H;
        _n_bands = _H({0, 0}, _p).rows();
    }

    void setParameters(const Parameters &p) { _p = p; }

    Parameters &params() { return _p; }

    Eigen::VectorXd eigenvals(const Point2D &k);
    Eigen::MatrixXcd eigenvecs(const Point2D &k);
    std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> eigen(const Point2D &k);

    void printBandStructure(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec);
    void printBandStructureSlice(const std::string &output_filename, const Eigen::VectorXd &k_vec, int axis, double k0 = 0.0);
    void printBerryCurvature(const std::string &output_filename, const H2D &dvxH, const H2D &dvyH, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec);
    void printGap(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec);
    void printGapAlongContour(const std::string &output_filename, const std::vector<Point2D> &contour);

    std::vector<std::vector<Point2D>> findFSContours(double E = 0.0, double dk = 1e-4, double eps = 1e-9, const Point2D kx_range = {-M_PI, 0.0}, std::size_t n_kx = 1001);

    Eigen::VectorXd calcChernNumbers(const Eigen::Vector2<std::size_t> &n, double kmax = 2.0 * M_PI);
    Eigen::VectorXd calcChernNumbersDenserCenter(std::size_t n_dense, std::size_t n_sparse, double k_val);
    Eigen::VectorXd calcChernNumbersWithCustomGrid(const Eigen::VectorXd &kx, const Eigen::VectorXd &ky);
    Eigen::VectorXd calcChernNumbersFromBC(const Eigen::Vector2<std::size_t> &n, const H2D &dvxH, const H2D &dvyH);

    H2D _H;
    Parameters _p;
    std::size_t _n_bands;

private:
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> _SAES;
};

#endif