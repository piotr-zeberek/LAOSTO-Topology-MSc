#include <iostream>
#include <iomanip>
#include <fstream>

#include "utils.h"
#include "System2D.h"

#include <unsupported/Eigen/KroneckerProduct>

// LAO/STO 001
int main2(int argc, char *argv[])
{
    Parameters p;

    p.tl = meV2au(875.0);
    p.th = meV2au(40.0);
    p.td = meV2au(40.0);
    p.delta_E = meV2au(47.0);
    p.delta_SO = meV2au(10.0);
    p.delta_RSO = meV2au(20.0);
    p.delta_SC = meV2au(0.2);
    p.g_Lande = 3;
    p.dx = nm2au(0.39);
    p.dy = nm2au(0.39);

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
        // account for the superconducting gap
        Eigen::MatrixXcd res = -p.delta_SC * Eigen::kroneckerProduct(sy, Eigen::kroneckerProduct(Eigen::Matrix3cd::Identity(), sy));

        const std::size_t Hk_dim = 6;

        res.topLeftCorner(Hk_dim, Hk_dim) = Hk(k, p);
        res.bottomRightCorner(Hk_dim, Hk_dim) = -Hk(-k, p).transpose();

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

    // complex delta testing
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(-3.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(5.0);

        // std::size_t n_k = 501;
        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -0.3, 0.3);
        // Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -0.3, 0.3);

        // sys.printAbsDelta("data/gap.dat", kx_vec, ky_vec);
    }

    // delta matrix for pairs
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(-3.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(5.0);

        // std::size_t n_k = 1001;
        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -0.5, 0.5);
        // Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -0.5, 0.5);

        // sys.printDeltaFromUnitaryTransformation("data/delta.dat", "data/DT.dat", kx_vec, ky_vec);
        // // sys.printAbsDelta("data/abs_delta.dat", kx_vec, ky_vec);

        // sys.printBandStructureSlice("data/BS.dat", kx_vec, 0, 0.0);

        // sys.setHamiltonian(Hk);
        // auto FS = sys.findFSContours(0.0, 1e-3, meV2au(1e-6));
        // sys.setHamiltonian(HBdG);

        // for (auto ic = 0; ic < FS.size(); ++ic)
        // {
        //     sys.printAbsDeltaAlongContour("data/gap_FS" + std::to_string(ic) + ".dat", FS[ic]);
        // }
    }

    // discrete ky
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(5.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(5.0);

        // std::size_t n_kx = 1001;
        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_kx, -0.8, 0.8);

        // std::size_t n_ky = 9;

        // sys.printBandStructureDiscreteky("data/BS.dat", kx_vec, n_ky);
        // sys.printBandStructureSlice("data/BS_slice.dat", kx_vec, 0, 0.0);
    }

    // abs delta vs Bz for given k
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0); // doesnt matter
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.0);

        // double Bz_max = T2au(5);
        // Eigen::VectorXd Bz_vec = Eigen::VectorXd::LinSpaced(10001, -Bz_max, Bz_max);

        // Point2D k{0.0, 0.0};

        // std::ofstream output_file("data/abs_delta_Bz.dat");

        // output_file << "# k = " << k.transpose() << std::endl;

        // for (auto Bz : Bz_vec)
        // {
        //     sys._p.Bz = Bz;
        //     output_file << Bz / T2au(1.0) << " " << sys.calcAbsDelta(k).transpose() / meV2au(1) << std::endl;
        // }
    }

    // BS
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(4.0);
        // sys._p.Bx = T2au(5.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.0);

        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(10001, -1.5, 1.5);

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

        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(10.0);

        // std::size_t n_dense = 10;  // number of k-points in each direction - even,
        // std::size_t n_sparse = 10; // number of k-points in each direction - even,

        // std::ofstream CN_file("data/CN.dat");

        // Eigen::VectorXd vec = Eigen::VectorXd::LinSpaced(101, -5.0, 5.0);

        // for (auto i = 0; i < vec.size(); ++i)
        // {
        //     sys._p.mu = meV2au(vec(i));
        //     CN_file << vec(i) << "\t" << sys.calcCNUsingWilsonLoop(n_dense, n_sparse, 1.0) << std::endl;
        //     std::cout << i + 1 << " out of " << vec.size() << std::endl;
        // }
    }

    // Chern numbers vs mu and B
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.0);

        // std::size_t n_dense = 20;  // number of k-points in each direction - even,
        // std::size_t n_sparse = 10; // number of k-points in each direction - even,


        // Eigen::VectorXd mu_vec = Eigen::VectorXd::LinSpaced(101, -5, 5);
        // Eigen::VectorXd Bz_vec = Eigen::VectorXd::LinSpaced(101, -5, 5);

        // std::ofstream mu_Bz_CN_file("data/mu_Bz_CN.dat");

        // for (auto i = 0; i < mu_vec.size(); ++i)
        // {
        //     for (auto j = 0; j < Bz_vec.size(); ++j)
        //     {
        //         sys._p.mu = meV2au(mu_vec(i));
        //         sys._p.Bz = T2au(Bz_vec(j));
        //         mu_Bz_CN_file << mu_vec(i) << "\t"
        //                       << Bz_vec(j) << "\t"
        //                       << sys.calcCNUsingWilsonLoop(n_dense, n_sparse, 0.6) << "\n";
        //         std::cout << i * Bz_vec.size() + j + 1 << " out of " << mu_vec.size() * Bz_vec.size() << std::endl;
        //     }
        // }
    }

    // Chern numbers vs mu and B - command line arguments version
    {
        // System2D sys(HBdG, p);

        // std::size_t n_dense = 5000;  // number of k-points in each direction - even,
        // std::size_t n_sparse = 1500; // number of k-points in each direction - even,

        // switch (argc)
        // {
        // case 2:
        //     sys._p.mu = meV2au(std::stod(argv[1]));
        //     break;
        // case 3:
        //     sys._p.mu = meV2au(std::stod(argv[1]));
        //     sys._p.Bz = T2au(std::stod(argv[2]));
        //     break;
        // case 5:
        //     sys._p.mu = meV2au(std::stod(argv[1]));
        //     sys._p.Bx = T2au(std::stod(argv[2]));
        //     sys._p.By = T2au(std::stod(argv[3]));
        //     sys._p.Bz = T2au(std::stod(argv[4]));
        //     break;
        // default:
        //     std::cerr << "Usage: " << argv[0] << " mu (meV)" << std::endl;
        //     std::cerr << "Usage: " << argv[0] << " mu (meV) Bz (T)" << std::endl;
        //     std::cerr << "Usage: " << argv[0] << " mu (meV) Bx By Bz (T)" << std::endl;
        //     return 1;
        // }

        // std::cout << "mu = " << sys._p.mu / meV2au(1.0) << " meV" << std::endl;
        // std::cout << "Bx = " << sys._p.Bx / T2au(1.0) << " T" << std::endl;
        // std::cout << "By = " << sys._p.By / T2au(1.0) << " T" << std::endl;
        // std::cout << "Bz = " << sys._p.Bz / T2au(1.0) << " T" << std::endl;

        // std::cout << "Chern numbers: " << sys.calcChernNumbersDenserCenter(n_dense, n_sparse, 0.5).transpose() << std::endl;
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

        // sys.printOrdinaryGap("data/gap.dat", kx_vec, kx_vec);
    }

    // FS contours and gap along them
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.1);

        // sys.setHamiltonian(Hk);
        // auto FS = sys.findFSContours();
        // sys.setHamiltonian(HBdG);

        // for (auto ic = 0; ic < FS.size(); ++ic)
        // {
        //     // sys.printOrdinaryGapAlongContour("data/gap_FS" + std::to_string(ic) + ".dat", FS[ic]);
        //     sys.printAbsDeltaAlongContour("data/gap_FS" + std::to_string(ic) + ".dat", FS[ic]);
        // }
    }

    return 0;
}

// check model
int main()
{
    Parameters p;

    p.m = 0.014;
    p.delta_SC = meV2au(0.2);
    p.g_Lande = -50;

    double dx = 8.0; // lattice constant in x in nm
    double dy = 8.0; // lattice constant in y in nm

    p.dx = nm2au(dx);
    p.dy = nm2au(dy);

    p.tx = 1.0 / (2.0 * p.m * p.dx * p.dx);
    p.ty = 1.0 / (2.0 * p.m * p.dy * p.dy);

    double alpha = 50.0; // Rashba coupling in meV nm
    p.delta_SO_x = meV2au(alpha / dx);
    p.delta_SO_y = meV2au(alpha / dy);

    auto Hk = [](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
    {
        // kinetic energy
        double ek = 2.0 * (p.tx * (1 - std::cos(k.x())) + p.ty * (1 - std::cos(k.y())));

        Eigen::MatrixXcd H0 = (ek - p.mu) * s0;

        // External magnetic field
        Eigen::MatrixXcd HB = 0.5 * p.g_Lande * 0.5 * (p.Bx * sx + p.By * sy + p.Bz * sz);

        // Rashba spin-orbit coupling
        Eigen::MatrixXcd HRSO = p.delta_SO_y * std::sin(k.y()) * sx - p.delta_SO_x * std::sin(k.x()) * sy;

        // final k-space hamiltonian
        return H0 + HB + HRSO;
    };

    auto HBdG = [Hk](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
    {
        // account for the superconducting gap
        Eigen::MatrixXcd res = -p.delta_SC * Eigen::kroneckerProduct(sy, sy);

        res.topLeftCorner(2, 2) = Hk(k, p);
        res.bottomRightCorner(2, 2) = -Hk(-k, p).transpose();

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

    // complex delta testing
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0);
        // sys._p.Bx = T2au(1.0);
        // sys._p.By = T2au(1.0);
        // sys._p.Bz = T2au(2.0);

        // // std::size_t n_k = 501;
        // // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -1.2, 1.2);
        // // Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -1.2, 1.2);

        // // sys.printGapFromUnitaryTransformation("data/gap.dat", kx_vec, ky_vec);

        // Point2D k{-0.1, 0};
        // Eigen::MatrixXcd H = sys._H(k, sys._p);
        // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> SAES(H);

        // std::cout << "E: " << SAES.eigenvalues().transpose() / meV2au(1) << std::endl;

        // Eigen::MatrixXcd U_tl = SAES.compute(sys._H(k, sys._p).topLeftCorner(2, 2)).eigenvectors();
        // std::cout << "Ep: " << SAES.eigenvalues().transpose() / meV2au(1) << std::endl;
        // Eigen::MatrixXcd U_br = SAES.compute(sys._H(k, sys._p).bottomRightCorner(2, 2)).eigenvectors().rowwise().reverse();
        // std::cout << "Em: " << -SAES.eigenvalues().transpose().reverse() / meV2au(1) << std::endl;

        // Eigen::MatrixXcd U = Eigen::MatrixXcd::Zero(4, 4);
        // U.topLeftCorner(2, 2) = U_tl;
        // U.bottomRightCorner(2, 2) = U_br;

        // H = (U.adjoint() * H * U).eval();

        // Eigen::MatrixXcd Ep = H.topLeftCorner(2, 2);
        // Eigen::MatrixXcd Em = -H.bottomRightCorner(2, 2);
        // Eigen::MatrixXcd D = H.topRightCorner(2, 2);
        // Eigen::MatrixXcd D_hc = H.bottomLeftCorner(2, 2);

        // std::ofstream output_file("data/det.dat");

        // double eps = meV2au(0.1);

        // for (double lambda = meV2au(-8); lambda <= meV2au(8); lambda += meV2au(0.0001))
        // {
        //     double detp = (-Em - lambda * s0 - D_hc * (Ep - lambda * s0).inverse() * D).determinant().real();
        //     double detm = (Ep - lambda * s0 - D * (-Em - lambda * s0).inverse() * D_hc).determinant().real();

        //     // bool flag = false;

        //     // for (int i = 0; i < Ep.rows(); ++i)
        //     // {
        //     //     if (std::abs(Ep(i, i) - lambda) < eps)
        //     //     {
        //     //         flag = true;
        //     //         break;
        //     //     }
        //     // }

        //     // if (flag)
        //     // {
        //     //     continue;
        //     // }

        //     output_file << lambda / meV2au(1) << " "
        //                 << detp / meV2au(1) << " "
        //                 << detm / meV2au(1) << std::endl;
        // }
    }

    // overlaps testing
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(1.0);
        // sys._p.Bx = T2au(4.0);
        // sys._p.By = T2au(1.0);
        // sys._p.Bz = T2au(0.5);

        // std::size_t n_k = 401;
        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -M_PI, M_PI);
        // Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -M_PI, M_PI);

        // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> SAES;

        // std::size_t b = 4;
        // std::size_t hb = b / 2;

        // Eigen::MatrixXd overlaps(hb, hb);

        // std::ofstream output_file("data/overlaps.dat");

        // for (auto kx : kx_vec)
        // {
        //     for (auto ky : ky_vec)
        //     {
        //         Point2D k{kx, ky};
        //         Eigen::MatrixXcd H = sys._H(k, sys._p);

        //         Eigen::MatrixXcd U_tl = SAES.compute(H.topLeftCorner(hb, hb)).eigenvectors();
        //         Eigen::MatrixXcd U_br = SAES.compute(H.bottomRightCorner(hb, hb)).eigenvectors();

        //         Eigen::MatrixXcd U = SAES.compute(H).eigenvectors();

        //         output_file << kx << " " << ky << " ";

        //         for (auto i = 0; i < b; ++i)
        //         {
        //             for (auto i_tl = 0; i_tl < hb; ++i_tl)
        //             {
        //                 for (auto i_br = 0; i_br < hb; ++i_br)
        //                 {
        //                     Eigen::VectorXcd v(b);
        //                     v << U_tl.col(i_tl), U_br.col(i_br);

        //                     overlaps(i_tl, i_br) = std::norm(U.col(i).dot(v));
        //                 }
        //             }

        //             output_file << overlaps.maxCoeff() << " ";
        //         }

        //         output_file << std::endl;
        //     }
        // }
    }

    // delta matrix with weird U
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.9);

        // std::size_t n_k = 401;
        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -1.2, 1.2);
        // Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -1.2, 1.2);

        // sys.printDeltaFromUnitaryTransformation("data/delta.dat", "data/DT.dat", kx_vec, ky_vec);
        // // sys.printAbsDelta("data/abs_delta.dat", kx_vec, ky_vec);

        // // sys.printBandStructureSlice("data/BS.dat", kx_vec, 0, 0.0);

        // sys.setHamiltonian(Hk);
        // auto FS = sys.findFSContours();
        // sys.setHamiltonian(HBdG);

        // for (auto ic = 0; ic < FS.size(); ++ic)
        // {
        //     sys.printOrdinaryGapAlongContour("data/gap_FS" + std::to_string(ic) + ".dat", FS[ic]);
        // }
    }

    // weird U
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.5);

        // std::size_t n_k = 401;
        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -M_PI, M_PI);
        // Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -M_PI, M_PI);

        // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> SAES;

        // std::size_t b = 4;
        // std::size_t hb = b / 2;

        // Point2D k{0, 0};
        // Eigen::MatrixXcd H = sys._H(k, sys._p);

        // Eigen::MatrixXcd C = kron(sx, s0);
        // Eigen::MatrixXcd D = Eigen::MatrixXcd::Zero(b, b);

        // SAES.compute(C);

        // for (auto i = 0; i < b; ++i)
        // {
        //     D(i, i) = SAES.eigenvalues()(i) > 0 ? 1 / std::sqrt(SAES.eigenvalues()(i)) : Eigen::dcomplex(0, 1 / std::sqrt(-SAES.eigenvalues()(i)));
        // }
        // Eigen::MatrixXcd U = SAES.eigenvectors() * D;
        // Eigen::MatrixXcd H_test = U.adjoint() * H * U;

        // // // PHS U
        // std::cout << C << std::endl;
        // std::cout << D << std::endl;
        // std::cout << U << std::endl;
        // std::cout << H_test / meV2au(1) << std::endl;
    }

    // discrete ky band structure
    {
        System2D sys(Hk, p);

        sys._p.mu = meV2au(0.0);
        sys._p.Bx = T2au(0.0);
        sys._p.By = T2au(0.0);
        sys._p.Bz = T2au(1.0);

        std::size_t n_kx = 5001;
        Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_kx, -1.0, 1.0);

        std::size_t n_ky = 21;

        sys.printBandStructureDiscreteky("data/BS.dat", kx_vec, n_ky);
        sys.printBandStructureSlice("data/BS_slice.dat", kx_vec, 0, 0.0);
    }

    // discrete hamiltonian energy vs parameters
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(4.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.0);

        // std::size_t n_kx = 250;
        // std::size_t n_ky = 9;
        // std::size_t n = 51;

        // double start = 0.0;
        // double end = 1.0;
        // Eigen::VectorXd vec = Eigen::VectorXd::LinSpaced(n, start, end);

        // std::ofstream output_file("data/energy.dat");

        // for (auto i = 0; i < n; ++i)
        // {
        //     sys._p.Bz = T2au(vec(i));
        //     // sys._p.mu = meV2au(vec(i));
        //     auto [evals, evecs] = sys.H_discrete_eigendecomposiiton(n_kx, n_ky);
        //     output_file << vec(i) << " " << evals.transpose() / meV2au(1) << std::endl;
        // }
    }

    // discrete system prob den
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.2);
        // sys._p.Bz = T2au(0.05);

        // std::size_t n_kx = 250;
        // std::size_t n_ky = 9;

        // double energy = meV2au(0.0);

        // sys.printProbDenDiscrete("data/prob_den.dat", n_kx, n_ky, energy);
    }

    // abs delta vs Bz for given k
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0); // doesnt matter
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.0);

        // double Bz_max = T2au(2);
        // Eigen::VectorXd Bz_vec = Eigen::VectorXd::LinSpaced(10001, -Bz_max, Bz_max);

        // Point2D k{0.0, 0.0};

        // std::ofstream output_file("data/abs_delta_Bz.dat");

        // output_file << "# k = " << k.transpose() << std::endl;

        // for (auto Bz : Bz_vec)
        // {
        //     sys._p.Bz = Bz;
        //     output_file << Bz / T2au(1.0) << " " << sys.calcAbsDelta(k).transpose() / meV2au(1) << std::endl;
        // }
    }

    // gap closing with changes in mu
    {
        // System2D sys(HBdG, p);

        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(5001, -M_PI, M_PI);

        // double mu_start = meV2au(-1.0);
        // double mu_end = meV2au(1.0);

        // Eigen::VectorXd mu_vec = Eigen::VectorXd::LinSpaced(9, mu_start, mu_end);
        // sys._p.Bz = T2au(0.4);

        // for (auto i = 0; i < mu_vec.size(); ++i)
        // {
        //     sys._p.mu = mu_vec(i);
        //     // sys.printBandStructureSlice("data/BS" + std::to_string(i) + ".dat", kx_vec, 0, 0.0);
        //     sys.printBandStructureDiscreteky("data/BS" + std::to_string(i) + ".dat", kx_vec, 5);
        // }
    }

    // Chern numbers vs parameter
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(4.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.0);

        // std::size_t n_dense = 10;  // number of k-points in each direction - even,
        // std::size_t n_sparse = 10; // number of k-points in each direction - even,

        // std::ofstream mu_Cs_file("data/CN.dat");

        // Eigen::VectorXd vec = Eigen::VectorXd::LinSpaced(101, -1.0, 1.0);

        // for (auto i = 0; i < vec.size(); ++i)
        // {
        //     sys._p.Bz = T2au(vec(i));
        //     // mu_Cs_file << mu_vec(i) / meV2au(1.0) << "\t" << sys.calcChernNumbersDenserCenter(n_dense, n_sparse, 0.5).transpose() << std::endl;
        //     // mu_Cs_file << mu_vec(i) / meV2au(1.0) << "\t" << sys.calcChernNumbersFromBC({1000, 1000}, dvxHBdG_num, dvyHBdG_num).transpose() << std::endl;
        //     mu_Cs_file << vec(i) << "\t" << sys.calcCNDiscreteky(n_dense, n_sparse, 0.8, 9) << std::endl;
        //     // mu_Cs_file << vec(i) << "\t" << sys.calcCNUsingWilsonLoop(n_dense, n_sparse, 0.5) << std::endl;
        //     std::cout << i + 1 << " out of " << vec.size() << std::endl;
        // }
    }

    // Chern numbers vs mu and Bz
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.0);

        // std::size_t n_dense = 20;  // number of k-points in each direction - even,
        // std::size_t n_sparse = 10; // number of k-points in each direction - even,

        // double mu_max = meV2au(1.0);
        // double Bz_max = T2au(0.5);

        // Eigen::VectorXd mu_vec = Eigen::VectorXd::LinSpaced(101, -2, 10);
        // Eigen::VectorXd Bz_vec = Eigen::VectorXd::LinSpaced(101, -1.0, 1.0);

        // std::ofstream mu_Bz_Cs_file("data/mu_Bz_CN.dat");
        // std::size_t n_ky = 21;

        // for (auto i = 0; i < mu_vec.size(); ++i)
        // {
        //     for (auto j = 0; j < Bz_vec.size(); ++j)
        //     {
        //         sys._p.mu = meV2au(mu_vec(i));
        //         sys._p.Bz = T2au(Bz_vec(j));
        //         // mu_Bz_Cs_file << mu_vec(i) << "\t"
        //         //               << Bz_vec(j) << "\t"
        //         //               << sys.calcCNUsingWilsonLoop(n_dense, n_sparse, 0.5) << "\n";
        //         mu_Bz_Cs_file << mu_vec(i) << "\t"
        //                       << Bz_vec(j) << "\t"
        //                       << sys.calcCNDiscreteky(n_dense, n_sparse, 1.0, n_ky) << "\n";
        //         std::cout << i * Bz_vec.size() + j + 1 << " out of " << mu_vec.size() * Bz_vec.size() << std::endl;
        //     }
        // }
    }

    // BC
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.3);

        // std::size_t n_k = 1001;
        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -0.5, 0.5);
        // Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -0.5, 0.5);

        // sys.printBerryCurvature("data/BC.dat", dvxHBdG_num, dvyHBdG_num, kx_vec, ky_vec);
        // sys.printBerryCurvature("data/BC_num2.dat", dvxHBdG_num2, dvyHBdG_num2, kx_vec, ky_vec);
    }

    // gap vs momentum
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.15);

        // std::size_t n_k = 501;
        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -0.3, 0.3);
        // Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -0.3, 0.3);

        // sys.printDeltaFromUnitaryTransformation("data/delta.dat", "data/DT.dat", kx_vec, kx_vec);
    }

    // FS contours and gap along them
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.3);

        // sys.setHamiltonian(Hk);
        // auto FS = sys.findFSContours();
        // sys.setHamiltonian(HBdG);

        // for (auto ic = 0; ic < FS.size(); ++ic)
        // {
        //     sys.printAbsDeltaAlongContour("data/gap_FS" + std::to_string(ic) + ".dat", FS[ic]);
        // }
    }

    return 0;
}