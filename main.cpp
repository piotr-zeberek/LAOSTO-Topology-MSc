#include <iostream>
#include <iomanip>
#include <fstream>

#include "utils.h"
#include "System2D.h"

#include <unsupported/Eigen/KroneckerProduct>

// LAO/STO 001
int main(int argc, char *argv[])
{
    Parameters p;

    p.tl = meV2au(875.0);
    p.th = meV2au(40.0);
    p.td = meV2au(40.0);
    p.delta_E = meV2au(47.0);
    p.delta_SO = meV2au(10.0);
    p.delta_RSO = meV2au(20.0);
    p.delta_SC = meV2au(0.2);
    p.mu = meV2au(0.0);
    p.g_Lande = 3;
    p.Bx = T2au(0.0);
    p.By = T2au(0.0);
    p.Bz = T2au(5.0);

    auto Hk = [](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
    {
        // kinetic energy
        double ek_xy = 2 * p.tl * (2 - std::cos(k.x()) - std::cos(k.y())) - p.delta_E;
        double ek_xz = 2 * p.tl * (1 - std::cos(k.x())) + 2 * p.th * (1 - std::cos(k.y()));
        double ek_yz = 2 * p.tl * (1 - std::cos(k.y())) + 2 * p.th * (1 - std::cos(k.x()));
        double ek_h = 2 * p.td * std::sin(k.x()) * std::sin(k.y());

        Eigen::Matrix3cd ek;
        ek << ek_xy, 0.0, 0.0,
            0.0, ek_xz, ek_h,
            0.0, ek_h, ek_yz;

        Eigen::MatrixXcd H0 = Eigen::kroneckerProduct(ek, s0);

        // atomic spin-orbit coupling
        Eigen::Matrix2cd z = Eigen::Matrix2cd::Zero();
        Eigen::MatrixXcd HSO(H0.rows(), H0.cols());

        HSO << z, 1i * sx, -1i * sy,
            -1i * sx, z, 1i * sz,
            1i * sy, -1i * sz, z;
        HSO *= p.delta_SO / 3.0;

        // Rashba spin-orbit coupling
        Eigen::Matrix3cd rso;
        rso << 0.0, 1i * std::sin(k.y()), 1i * std::sin(k.x()),
            -1i * std::sin(k.y()), 0.0, 0.0,
            -1i * std::sin(k.x()), 0.0, 0.0;
        Eigen::MatrixXcd HRSO = p.delta_RSO * Eigen::kroneckerProduct(rso, s0);

        // External magnetic field
        Eigen::Matrix3cd I3 = Eigen::Matrix3cd::Identity();

        Eigen::MatrixXcd HBx = p.Bx * (Eigen::kroneckerProduct(Lx, s0) + 0.5 * p.g_Lande * Eigen::kroneckerProduct(I3, sx));
        Eigen::MatrixXcd HBy = p.By * (Eigen::kroneckerProduct(Ly, s0) + 0.5 * p.g_Lande * Eigen::kroneckerProduct(I3, sy));
        Eigen::MatrixXcd HBz = p.Bz * (Eigen::kroneckerProduct(Lz, s0) + 0.5 * p.g_Lande * Eigen::kroneckerProduct(I3, sz));

        // final k-space hamiltonian
        return H0 + HSO + HRSO + 0.5 * (HBx + HBy + HBz) - p.mu * Eigen::MatrixXcd::Identity(H0.rows(), H0.cols());
    };

    auto HBdG = [Hk](const Eigen::Vector2d &k, const Parameters &p) -> Eigen::MatrixXcd
    {
        Eigen::MatrixXcd upper_left = Hk(k, p);
        Eigen::MatrixXcd lower_right = -Hk(-k, p).transpose();

        // account for the superconducting gap
        Eigen::MatrixXcd res = -p.delta_SC * Eigen::kroneckerProduct(sy, Eigen::kroneckerProduct(Eigen::Matrix3cd::Identity(), sy));

        res.topLeftCorner(upper_left.rows(), upper_left.cols()) = upper_left;
        res.bottomRightCorner(lower_right.rows(), lower_right.cols()) = lower_right;

        return res;
    };

    double dk = 1e-4;

    auto dvxHBdG_num = [HBdG, dk](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
    { return (HBdG(k + Point2D{dk, 0.0}, p) - HBdG(k - Point2D{dk, 0.0}, p)) / (2.0 * dk); };

    auto dvyHBdG_num = [HBdG, dk](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
    { return (HBdG(k + Point2D{0.0, dk}, p) - HBdG(k - Point2D{0.0, dk}, p)) / (2.0 * dk); };

    auto dvxHBdG_num2 = [HBdG, dk](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
    { return (-HBdG(k + Point2D{2.0 * dk, 0.0}, p) + 8.0 * HBdG(k + Point2D{dk, 0.0}, p) - 8.0 * HBdG(k - Point2D{dk, 0.0}, p) + HBdG(k - Point2D{2.0 * dk, 0.0}, p)) / (12.0 * dk); };

    auto dvyHBdG_num2 = [HBdG, dk](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
    { return (-HBdG(k + Point2D{0.0, 2.0 * dk}, p) + 8.0 * HBdG(k + Point2D{0.0, dk}, p) - 8.0 * HBdG(k - Point2D{0.0, dk}, p) + HBdG(k - Point2D{0.0, 2.0 * dk}, p)) / (12.0 * dk); };

    // BS
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(-2.0);
        // sys._p.Bz = T2au(2.0);

        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(10001, -0.5, 0.5);

        // sys.printBandStructureSlice("data/BS.dat", kx_vec, 0, 0.0); // slice
        // // sys.printBandStructure("data/BS.dat", kx_vec, kx_vec); // 2D
    }

    // gap closing with changes in mu
    {
        // System2D sys(HBdG, p);

        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(10001, -0.5, 0.5);

        // double mu_start = meV2au(-2.5);
        // double mu_end = meV2au(-2.1);

        // Eigen::VectorXd mu_vec = Eigen::VectorXd::LinSpaced(9, mu_start, mu_end);

        // for (auto i = 0; i < mu_vec.size(); ++i)
        // {
        //     sys._p.mu = mu_vec(i);
        //     sys.printBandStructureSlice("data/BS" + std::to_string(i) + ".dat", kx_vec, 0, 0.0);
        // }
    }

    // Chern numbers vs mu
    {
        // System2D sys(HBdG, p);

        // std::size_t n_dense = 4000; // number of k-points in each direction - even,
        // std::size_t n_sparse = 1000; // number of k-points in each direction - even,

        // std::ofstream mu_Cs_file("data/mu_Cs.dat");

        // double mu_max = meV2au(4.0);
        // Eigen::VectorXd mu_vec = Eigen::VectorXd::LinSpaced(9, -mu_max, mu_max);

        // for (auto i = 0; i < mu_vec.size(); ++i)
        // {
        //     sys._p.mu = mu_vec(i);
        //     mu_Cs_file << mu_vec(i) / meV2au(1.0) << "\t" << sys.calcChernNumbersDenserCenter(n_dense, n_sparse, 0.5).transpose() << std::endl;
        //     std::cout << i + 1 << " out of " << mu_vec.size() << std::endl;
        // }
    }

    // Chern numbers vs mu and B
    {
        // System2D sys(HBdG, p);

        // std::size_t n_dense = 4000; // number of k-points in each direction - even,
        // std::size_t n_sparse = 1000; // number of k-points in each direction - even,

        // double mu_max = meV2au(4.0);
        // double Bz_max = T2au(5);

        // Eigen::VectorXd mu_vec = Eigen::VectorXd::LinSpaced(5, -mu_max, mu_max);
        // Eigen::VectorXd Bz_vec = Eigen::VectorXd::LinSpaced(6, -Bz_max, Bz_max);

        // std::ofstream mu_Bz_Cs_file("data/mu_Bz_Cs.dat");

        // for (auto i = 0; i < mu_vec.size(); ++i)
        // {
        //     for (auto j = 0; j < Bz_vec.size(); ++j)
        //     {
        //         sys._p.mu = mu_vec(i);
        //         sys._p.Bz = Bz_vec(j);
        //         mu_Bz_Cs_file << mu_vec(i) / meV2au(1.0) << "\t"
        //                       << Bz_vec(j) / T2au(1.0) << "\t"
        //                       << sys.calcChernNumbersDenserCenter(n_dense, n_sparse, 0.45).transpose() << std::endl;
        //         std::cout << i * Bz_vec.size() + j + 1 << " out of " << mu_vec.size() * Bz_vec.size() << std::endl;
        //     }
        // }
    }

    // Chern numbers vs mu and B - command line arguments version
    {
        System2D sys(HBdG, p);

        std::size_t n_dense = 5000;  // number of k-points in each direction - even,
        std::size_t n_sparse = 1500; // number of k-points in each direction - even,

        switch (argc)
        {
        case 2:
            sys._p.mu = meV2au(std::stod(argv[1]));
            break;
        case 3:
            sys._p.mu = meV2au(std::stod(argv[1]));
            sys._p.Bz = T2au(std::stod(argv[2]));
            break;
        case 5:
            sys._p.mu = meV2au(std::stod(argv[1]));
            sys._p.Bx = T2au(std::stod(argv[2]));
            sys._p.By = T2au(std::stod(argv[3]));
            sys._p.Bz = T2au(std::stod(argv[4]));
            break;
        default:
            std::cerr << "Usage: " << argv[0] << " mu (meV)" << std::endl;
            std::cerr << "Usage: " << argv[0] << " mu (meV) Bz (T)" << std::endl;
            std::cerr << "Usage: " << argv[0] << " mu (meV) Bx By Bz (T)" << std::endl;
            return 1;
        }

        std::cout << "mu = " << sys._p.mu / meV2au(1.0) << " meV" << std::endl;
        std::cout << "Bx = " << sys._p.Bx / T2au(1.0) << " T" << std::endl;
        std::cout << "By = " << sys._p.By / T2au(1.0) << " T" << std::endl;
        std::cout << "Bz = " << sys._p.Bz / T2au(1.0) << " T" << std::endl;

        std::cout << "Chern numbers: " << sys.calcChernNumbersDenserCenter(n_dense, n_sparse, 0.5).transpose() << std::endl;
    }

    // BC
    {
        // System2D sys(HBdG, p);
        // sys._p.mu = meV2au(-2);
        // sys._p.Bz = T2au(5);

        // std::size_t n_k = 1001;
        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -0.5, 0.5);
        // Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -0.5, 0.5);

        // sys.printBerryCurvature("data/BC.dat", dvxHBdG_num, dvyHBdG_num, kx_vec, ky_vec);
        // // sys.printBerryCurvature("data/BC_num2.dat", dvxHBdG_num2, dvyHBdG_num2, kx_vec, ky_vec);
    }

    // gap vs momentum
    {
        // System2D sys(HBdG, p);

        // std::size_t n_k = 1001;
        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -1.5, 1.5);
        // Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -1.5, 1.5);

        // sys._p.mu = meV2au(0.0);

        // sys.printGap("data/gap.dat", kx_vec, kx_vec);
    }

    // FS contours and gap along them
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0);
        // sys.setHamiltonian(Hk);
        // auto FS = sys.findFSContours();
        // sys.setHamiltonian(HBdG);

        // for (auto ic = 0; ic < FS.size(); ++ic)
        // {
        //     sys.printGapAlongContour("data/gap_FS" + std::to_string(ic) + ".dat", FS[ic]);
        // }
    }

    return 0;
}

// check model
int main2()
{
    Parameters p;

    p.m = 0.014;
    p.alpha = meV2au(50.0);
    p.delta_SC = meV2au(0.2);
    p.g_Lande = -50;
    p.Bx = T2au(0.0);
    p.By = T2au(0.0);
    p.Bz = T2au(2.0);

    p.a = nm2au(8.0);

    p.t = 1.0 / (2.0 * p.m * p.a * p.a);
    p.mu = meV2au(0.0);

    // std::cout << "t = " << p.t / meV2au(1.0) << " meV" << std::endl;

    auto Hk = [](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
    {
        // kinetic energy
        double ek = 4 * p.t - 2.0 * p.t * (std::cos(k.x()) + std::cos(k.y()));

        Eigen::MatrixXcd H0 = ek * s0;

        // Rashba spin-orbit coupling
        Eigen::MatrixXcd HRSO = p.alpha * (std::sin(k.y()) * sx - std::sin(k.x()) * sy);

        // External magnetic field
        Eigen::MatrixXcd HBx = p.Bx * sx;
        Eigen::MatrixXcd HBy = p.By * sy;
        Eigen::MatrixXcd HBz = p.Bz * sz;

        // final k-space hamiltonian
        return H0 + HRSO + 0.5 * p.g_Lande * 0.5 * (HBx + HBy + HBz) - p.mu * s0; // mu_Bohr = 0.5 in au
    };

    auto HBdG = [Hk](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
    {
        Eigen::MatrixXcd upper_left = Hk(k, p);
        Eigen::MatrixXcd lower_right = -Hk(-k, p).transpose();

        // account for the superconducting gap
        Eigen::MatrixXcd res = -p.delta_SC * Eigen::kroneckerProduct(sy, sy);

        res.topLeftCorner(upper_left.rows(), upper_left.cols()) = upper_left;
        res.bottomRightCorner(lower_right.rows(), lower_right.cols()) = lower_right;

        return res;
    };

    // numeryczne pochodne
    double dk = 1e-4;

    auto dvxHBdG_num = [HBdG, dk](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
    { return (HBdG(k + Point2D{dk, 0.0}, p) - HBdG(k - Point2D{dk, 0.0}, p)) / (2.0 * dk); };

    auto dvyHBdG_num = [HBdG, dk](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
    { return (HBdG(k + Point2D{0.0, dk}, p) - HBdG(k - Point2D{0.0, dk}, p)) / (2.0 * dk); };

    auto dvxHBdG_num2 = [HBdG, dk](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
    { return (-HBdG(k + Point2D{2.0 * dk, 0.0}, p) + 8.0 * HBdG(k + Point2D{dk, 0.0}, p) - 8.0 * HBdG(k - Point2D{dk, 0.0}, p) + HBdG(k - Point2D{2.0 * dk, 0.0}, p)) / (12.0 * dk); };

    auto dvyHBdG_num2 = [HBdG, dk](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
    { return (-HBdG(k + Point2D{0.0, 2.0 * dk}, p) + 8.0 * HBdG(k + Point2D{0.0, dk}, p) - 8.0 * HBdG(k - Point2D{0.0, dk}, p) + HBdG(k - Point2D{0.0, 2.0 * dk}, p)) / (12.0 * dk); };

    // BS
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(-1.0);

        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(2001, -M_PI, M_PI);

        // sys.printBandStructureSlice("data/BS.dat", kx_vec, 0, 0.0); // slice
        // // sys.printBandStructure("data/BS.dat", kx_vec, kx_vec); // 2D
    }

    // gap closing with changes in mu
    {
        // System2D sys(HBdG, p);

        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(10001, -1.5, 1.5);

        // double mu_start = meV2au(-1.0);
        // double mu_end = meV2au(1.0);

        // Eigen::VectorXd mu_vec = Eigen::VectorXd::LinSpaced(9, mu_start, mu_end);

        // for (auto i = 0; i < mu_vec.size(); ++i)
        // {
        //     sys._p.mu = mu_vec(i);
        //     sys.printBandStructureSlice("data/BS" + std::to_string(i) + ".dat", kx_vec, 0, 0.0);
        // }
    }

    // Chern numbers vs mu
    {
        // System2D sys(HBdG, p);

        // std::size_t n_dense = 1000; // number of k-points in each direction - even,
        // std::size_t n_sparse = 500; // number of k-points in each direction - even,

        // std::ofstream mu_Cs_file("data/mu_Cs.dat");

        // double mu_max = meV2au(4.0);
        // Eigen::VectorXd mu_vec = Eigen::VectorXd::LinSpaced(201, -mu_max, mu_max);

        // for (auto i = 0; i < mu_vec.size(); ++i)
        // {
        //     sys._p.mu = mu_vec(i);
        //     mu_Cs_file << mu_vec(i) / meV2au(1.0) << "\t" << sys.calcChernNumbersDenserCenter(n_dense, n_sparse, 1.2).transpose() << std::endl;
        //     std::cout << i + 1 << " out of " << mu_vec.size() << std::endl;
        // }
    }

    // Chern numbers vs mu and Bz
    {
        // System2D sys(HBdG, p);

        // std::size_t n_dense = 1000; // number of k-points in each direction - even,
        // std::size_t n_sparse = 500; // number of k-points in each direction - even,

        // double mu_max = meV2au(1.0);
        // double Bz_max = T2au(0.5);

        // Eigen::VectorXd mu_vec = Eigen::VectorXd::LinSpaced(201, -mu_max, mu_max);
        // Eigen::VectorXd Bz_vec = Eigen::VectorXd::LinSpaced(201, -Bz_max, Bz_max);

        // std::ofstream mu_Bz_Cs_file("data/mu_Bz_Cs.dat");

        // for (auto i = 0; i < mu_vec.size(); ++i)
        // {
        //     for (auto j = 0; j < Bz_vec.size(); ++j)
        //     {
        //         sys._p.mu = mu_vec(i);
        //         sys._p.Bz = Bz_vec(j);
        //         mu_Bz_Cs_file << mu_vec(i) / meV2au(1.0) << "\t"
        //                       << Bz_vec(j) / T2au(1.0) << "\t"
        //                       << sys.calcChernNumbersDenserCenter(n_dense, n_sparse, 1.2).transpose() << "\n";
        //         std::cout << i * Bz_vec.size() + j + 1 << " out of " << mu_vec.size() * Bz_vec.size() << std::endl;
        //     }
        // }
    }

    // BC
    {
        // System2D sys(HBdG, p);
        // sys._p.mu = meV2au(-2.89);

        // std::size_t n_k = 2001;
        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -M_PI, M_PI);
        // Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -M_PI, M_PI);

        // // analityczne pochodne - jest zle chyba
        // sys.printBerryCurvature("data/BC.dat", dvxHBdG, dvyHBdG, kx_vec, ky_vec);

        // sys.printBerryCurvature("data/BC_num.dat", dvxHBdG_num, dvyHBdG_num, kx_vec, ky_vec);
        // sys.printBerryCurvature("data/BC_num2.dat", dvxHBdG_num2, dvyHBdG_num2, kx_vec, ky_vec);
    }

    // gap vs momentum
    {
        // System2D sys(HBdG, p);

        // std::size_t n_k = 1001;
        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -1.5, 1.5);
        // Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -1.5, 1.5);

        // sys._p.mu = meV2au(0.0);

        // sys.printGap("data/gap.dat", kx_vec, kx_vec);
    }

    // FS contours and gap along them
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0);
        // sys.setHamiltonian(Hk);
        // auto FS = sys.findFSContours();
        // sys.setHamiltonian(HBdG);

        // for (auto ic = 0; ic < FS.size(); ++ic)
        // {
        //     sys.printGapAlongContour("data/gap_FS" + std::to_string(ic) + ".dat", FS[ic]);
        // }
    }

    return 0;
}