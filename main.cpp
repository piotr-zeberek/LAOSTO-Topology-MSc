#include <iostream>
#include <iomanip>
#include <fstream>

#include "ToyModelElementBased.h"
#include "ToyModelMatrixBased.h"

#include "LAOSTOMatrixBased.h"

#include "System2DCalculations.h"   
#include "System2DCalculationsPrinter.h"

// // LAO/STO 001
// int main2(int argc, char *argv[])
// {
//     Parameters p;

//     p.tl = meV2au(875.0);
//     p.th = meV2au(40.0);
//     p.td = meV2au(40.0);
//     p.delta_E = meV2au(47.0);
//     p.delta_SO = meV2au(10.0);
//     p.delta_RSO = meV2au(20.0);
//     p.delta_SC = meV2au(0.2);
//     p.g_Lande = 3;
//     p.dx = nm2au(0.39);
//     p.dy = nm2au(0.39);

//     auto Hk = [](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
//     {
//         // kinetic energy
//         double ek_xy = 2 * p.tl * (2 - std::cos(k.x()) - std::cos(k.y())) - p.delta_E;
//         double ek_xz = 2 * p.tl * (1 - std::cos(k.x())) + 2 * p.th * (1 - std::cos(k.y()));
//         double ek_yz = 2 * p.tl * (1 - std::cos(k.y())) + 2 * p.th * (1 - std::cos(k.x()));
//         double ek_h = 2 * p.td * std::sin(k.x()) * std::sin(k.y());

//         Eigen::Matrix3cd ek;
//         ek << ek_xy, 0.0, 0.0,
//             0.0, ek_xz, ek_h,
//             0.0, ek_h, ek_yz;

//         Eigen::MatrixXcd H0 = Eigen::kroneckerProduct(ek, s0);

//         // atomic spin-orbit coupling
//         Eigen::Matrix2cd z = Eigen::Matrix2cd::Zero();
//         Eigen::MatrixXcd HSO(H0.rows(), H0.cols());

//         HSO << z, 1i * sx, -1i * sy,
//             -1i * sx, z, 1i * sz,
//             1i * sy, -1i * sz, z;
//         HSO *= p.delta_SO / 3.0;

//         // Rashba spin-orbit coupling
//         Eigen::Matrix3cd rso;
//         rso << 0.0, 1i * std::sin(k.y()), 1i * std::sin(k.x()),
//             -1i * std::sin(k.y()), 0.0, 0.0,
//             -1i * std::sin(k.x()), 0.0, 0.0;
//         Eigen::MatrixXcd HRSO = p.delta_RSO * Eigen::kroneckerProduct(rso, s0);

//         // External magnetic field
//         Eigen::Matrix3cd I3 = Eigen::Matrix3cd::Identity();

//         Eigen::MatrixXcd HBx = p.Bx * (Eigen::kroneckerProduct(Lx, s0) + 0.5 * p.g_Lande * Eigen::kroneckerProduct(I3, sx));
//         Eigen::MatrixXcd HBy = p.By * (Eigen::kroneckerProduct(Ly, s0) + 0.5 * p.g_Lande * Eigen::kroneckerProduct(I3, sy));
//         Eigen::MatrixXcd HBz = p.Bz * (Eigen::kroneckerProduct(Lz, s0) + 0.5 * p.g_Lande * Eigen::kroneckerProduct(I3, sz));

//         // final k-space hamiltonian
//         return H0 + HSO + HRSO + 0.5 * (HBx + HBy + HBz) - p.mu * Eigen::MatrixXcd::Identity(H0.rows(), H0.cols());
//     };

//     auto HBdG = [Hk](const Eigen::Vector2d &k, const Parameters &p) -> Eigen::MatrixXcd
//     {
//         // account for the superconducting gap
//         Eigen::MatrixXcd res = -p.delta_SC * Eigen::kroneckerProduct(sy, Eigen::kroneckerProduct(Eigen::Matrix3cd::Identity(), sy));

//         const std::size_t Hk_dim = 6;

//         res.topLeftCorner(Hk_dim, Hk_dim) = Hk(k, p);
//         res.bottomRightCorner(Hk_dim, Hk_dim) = -Hk(-k, p).transpose();

//         return res;
//     };

//     double dk = 1e-4;

//     auto dvxHBdG_num = [HBdG, dk](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
//     { return (HBdG(k + Point2D{dk, 0.0}, p) - HBdG(k - Point2D{dk, 0.0}, p)) / (2.0 * dk); };

//     auto dvyHBdG_num = [HBdG, dk](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
//     { return (HBdG(k + Point2D{0.0, dk}, p) - HBdG(k - Point2D{0.0, dk}, p)) / (2.0 * dk); };

//     auto dvxHBdG_num2 = [HBdG, dk](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
//     { return (-HBdG(k + Point2D{2.0 * dk, 0.0}, p) + 8.0 * HBdG(k + Point2D{dk, 0.0}, p) - 8.0 * HBdG(k - Point2D{dk, 0.0}, p) + HBdG(k - Point2D{2.0 * dk, 0.0}, p)) / (12.0 * dk); };

//     auto dvyHBdG_num2 = [HBdG, dk](const Point2D &k, const Parameters &p) -> Eigen::MatrixXcd
//     { return (-HBdG(k + Point2D{0.0, 2.0 * dk}, p) + 8.0 * HBdG(k + Point2D{0.0, dk}, p) - 8.0 * HBdG(k - Point2D{0.0, dk}, p) + HBdG(k - Point2D{0.0, 2.0 * dk}, p)) / (12.0 * dk); };

//     // // delta vs momentum
//     // {
//     //     System2D sys(HBdG, p);

//     //     sys._p.mu = meV2au(0.0);
//     //     sys._p.Bx = T2au(0.0);
//     //     sys._p.By = T2au(0.0);
//     //     sys._p.Bz = T2au(0.05);

//     //     std::size_t n_k = 201;
//     //     double k_max = 0.5;

//     //     Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -k_max, k_max);
//     //     Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -k_max, k_max);

//     //     sys.printAbsDelta("data/abs_delta.dat", kx_vec, ky_vec);
//     //     sys.printDeltaFromUnitaryTransformation("data/delta.dat", "data/DT.dat", kx_vec, ky_vec);
//     // }

//     // discrete ky
//     {
//         // System2D sys(HBdG, p);

//         // sys._p.mu = meV2au(0.0);
//         // sys._p.Bx = T2au(0.0);
//         // sys._p.By = T2au(0.0);
//         // sys._p.Bz = T2au(0.1);

//         // std::size_t n_kx = 1001;
//         // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_kx, -0.8, 0.8);

//         // std::size_t n_ky = 9;

//         // sys.printBandStructureDiscreteky("data/BS.dat", kx_vec, n_ky);
//         // sys.printBandStructureSlice("data/BS_slice.dat", kx_vec, 0, 0.0);
//     }

//     // abs delta vs Bz for given k
//     {
//         // System2D sys(HBdG, p);

//         // sys._p.mu = meV2au(0.0); // doesnt matter
//         // sys._p.Bx = T2au(0.0);
//         // sys._p.By = T2au(0.0);
//         // sys._p.Bz = T2au(0.0);

//         // double Bz_max = T2au(5);
//         // Eigen::VectorXd Bz_vec = Eigen::VectorXd::LinSpaced(10001, -Bz_max, Bz_max);

//         // Point2D k{0.0, 0.0};

//         // std::ofstream output_file("data/abs_delta_Bz.dat");

//         // output_file << "# k = " << k.transpose() << std::endl;

//         // for (auto Bz : Bz_vec)
//         // {
//         //     sys._p.Bz = Bz;
//         //     output_file << Bz / T2au(1.0) << " " << sys.calcAbsDelta(k).transpose() / meV2au(1) << std::endl;
//         // }
//     }

//     // BS
//     {
//         // System2D sys(HBdG, p);

//         // sys._p.mu = meV2au(4.0);
//         // sys._p.Bx = T2au(5.0);
//         // sys._p.By = T2au(0.0);
//         // sys._p.Bz = T2au(0.0);

//         // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(10001, -1.5, 1.5);

//         // sys.printBandStructureSlice("data/BS.dat", kx_vec, 0, 0.0); // slice
//         // // sys.printBandStructure("data/BS.dat", kx_vec, kx_vec); // 2D
//     }

//     // gap closing with changes in mu
//     {
//         // System2D sys(HBdG, p);

//         // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(10001, -0.5, 0.5);

//         // double mu_start = meV2au(-2.5);
//         // double mu_end = meV2au(-2.1);

//         // Eigen::VectorXd mu_vec = Eigen::VectorXd::LinSpaced(9, mu_start, mu_end);

//         // for (auto i = 0; i < mu_vec.size(); ++i)
//         // {
//         //     sys._p.mu = mu_vec(i);
//         //     sys.printBandStructureSlice("data/BS" + std::to_string(i) + ".dat", kx_vec, 0, 0.0);
//         // }
//     }

//     // Chern numbers vs mu
//     {
//         // System2D sys(HBdG, p);

//         // sys._p.Bx = T2au(0.0);
//         // sys._p.By = T2au(0.0);
//         // sys._p.Bz = T2au(5.0);

//         // std::size_t n_dense = 500;  // number of k-points in each direction - even,
//         // std::size_t n_sparse = 100; // number of k-points in each direction - even,

//         // std::ofstream CN_file("data/CN.dat");

//         // Eigen::VectorXd vec = Eigen::VectorXd::LinSpaced(21, -4.0, -2.0);

//         // for (auto i = 0; i < vec.size(); ++i)
//         // {
//         //     sys._p.mu = meV2au(vec(i));
//         //     CN_file << vec(i)
//         //             // << "\t" << sys.calcCNUsingWilsonLoop(n_dense, n_sparse, 1.0) << std::endl;
//         //             << "\t" << sys.calcChernNumberUsingMatrixBerryCurvatureTrace(n_dense, n_sparse, 1.0, 1000) << std::endl;
//         //     std::cout << i + 1 << " out of " << vec.size() << std::endl;
//         // }
//     }

//     // Chern numbers vs mu and B
//     {
//         // System2D sys(HBdG, p);

//         // sys._p.mu = meV2au(0.0);
//         // sys._p.Bx = T2au(0.0);
//         // sys._p.By = T2au(0.0);
//         // sys._p.Bz = T2au(0.5);

//         // std::size_t n_dense = 50; // number of k-points in each direction - even,
//         // std::size_t n_sparse = 8; // number of k-points in each direction - even,

//         // double mu_max = meV2au(1.0);
//         // double Bz_max = T2au(0.5);

//         // Eigen::VectorXd mu_vec = Eigen::VectorXd::LinSpaced(21, -1, 1);
//         // Eigen::VectorXd B_vec = Eigen::VectorXd::LinSpaced(21, -0.5, 0.5);

//         // std::ofstream mu_B_CN_file("data/mu_B_CN.dat");
//         // std::size_t n_ky = 21;

//         // for (auto i = 0; i < mu_vec.size(); ++i)
//         // {
//         //     for (auto j = 0; j < B_vec.size(); ++j)
//         //     {
//         //         sys._p.mu = meV2au(mu_vec(i));
//         //         sys._p.By = T2au(B_vec(j));
//         //         mu_B_CN_file << mu_vec(i) << "\t"
//         //                      << B_vec(j) << "\t"
//         //                      //   << sys.calcCNDiscreteky(n_dense, n_sparse, 1.0, n_ky) << "\n";
//         //                      //   << sys.calcCNUsingWilsonLoop(n_dense, n_sparse, 0.5) << std::endl;
//         //                      << sys.calcChernNumberUsingMatrixBerryCurvatureTrace(n_dense, n_sparse, 0.5) << std::endl;
//         //         std::cout << i * B_vec.size() + j + 1 << " out of " << mu_vec.size() * B_vec.size() << std::endl;
//         //     }
//         // }
//     }

//     // Chern numbers vs mu and B - command line arguments version
//     {
//         // System2D sys(HBdG, p);

//         // std::size_t n_dense = 5000;  // number of k-points in each direction - even,
//         // std::size_t n_sparse = 1500; // number of k-points in each direction - even,

//         // switch (argc)
//         // {
//         // case 2:
//         //     sys._p.mu = meV2au(std::stod(argv[1]));
//         //     break;
//         // case 3:
//         //     sys._p.mu = meV2au(std::stod(argv[1]));
//         //     sys._p.Bz = T2au(std::stod(argv[2]));
//         //     break;
//         // case 5:
//         //     sys._p.mu = meV2au(std::stod(argv[1]));
//         //     sys._p.Bx = T2au(std::stod(argv[2]));
//         //     sys._p.By = T2au(std::stod(argv[3]));
//         //     sys._p.Bz = T2au(std::stod(argv[4]));
//         //     break;
//         // default:
//         //     std::cerr << "Usage: " << argv[0] << " mu (meV)" << std::endl;
//         //     std::cerr << "Usage: " << argv[0] << " mu (meV) Bz (T)" << std::endl;
//         //     std::cerr << "Usage: " << argv[0] << " mu (meV) Bx By Bz (T)" << std::endl;
//         //     return 1;
//         // }

//         // std::cout << "mu = " << sys._p.mu / meV2au(1.0) << " meV" << std::endl;
//         // std::cout << "Bx = " << sys._p.Bx / T2au(1.0) << " T" << std::endl;
//         // std::cout << "By = " << sys._p.By / T2au(1.0) << " T" << std::endl;
//         // std::cout << "Bz = " << sys._p.Bz / T2au(1.0) << " T" << std::endl;

//         // std::cout << "Chern numbers: " << sys.calcChernNumbersDenserCenter(n_dense, n_sparse, 0.5).transpose() << std::endl;
//     }

//     // BC
//     {
//         // System2D sys(HBdG, p);

//         // sys._p.mu = meV2au(3.2);
//         // sys._p.Bx = T2au(0.0);
//         // sys._p.By = T2au(0.0);
//         // sys._p.Bz = T2au(5);

//         // std::size_t n_k = 201;
//         // double k_max = 0.5;

//         // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -k_max, k_max);
//         // Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -k_max, k_max);

//         // sys.printBerryCurvature("data/BC.dat", kx_vec, ky_vec);
//         // sys.printBerryCurvatureFromWilsonLoop("data/BC_WL.dat", kx_vec, ky_vec);
//     }

//     // wilson loop spectrum
//     {
//         // System2D sys(HBdG, p);

//         // sys._p.mu = meV2au(-3.0);
//         // sys._p.Bx = T2au(0.0);
//         // sys._p.By = T2au(0.0);
//         // sys._p.Bz = T2au(5);

//         // sys.printWilsonLoopSpectrum("data/WL.dat", 201, 50001, 0.3);
//     }

//     // delta vs momentum
//     {
//         System2D sys(HBdG, p);

//         sys._p.mu = meV2au(-3.0);
//         sys._p.Bx = T2au(0.0);
//         sys._p.By = T2au(0.0);
//         sys._p.Bz = T2au(5);

//         std::size_t n_k = 501;
//         Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -0.5, 0.5);
//         Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -0.5, 0.5);

//         sys.printAbsDelta("data/gap.dat", kx_vec, kx_vec);
//         sys.printDeltaFromUnitaryTransformation("data/delta.dat", "data/DT.dat", kx_vec, kx_vec);
//     }

//     // FS contours and gap along them
//     {
//         // System2D sys(HBdG, p);

//         // sys._p.mu = meV2au(0.0);
//         // sys._p.Bx = T2au(0.0);
//         // sys._p.By = T2au(0.0);
//         // sys._p.Bz = T2au(0.1);

//         // sys.setHamiltonian(Hk);
//         // auto FS = sys.findFSContours();
//         // sys.setHamiltonian(HBdG);

//         // for (auto ic = 0; ic < FS.size(); ++ic)
//         // {
//         //     // sys.printOrdinaryGapAlongContour("data/gap_FS" + std::to_string(ic) + ".dat", FS[ic]);
//         //     sys.printAbsDeltaAlongContour("data/gap_FS" + std::to_string(ic) + ".dat", FS[ic]);
//         // }
//     }

//     return 0;
// }

// check model
int main()
{
    ToyModelMatrixBased sys;
    System2DCalculations calc(sys);
    System2DCalculationsPrinter printer(calc);

    ToyModelElementBased sys2;
    System2DCalculations calc2(sys2);
    System2DCalculationsPrinter printer2(calc2);

    // LAOSTO sys;
    // System2DCalculations calc(sys);
    // System2DCalculationsPrinter printer(calc);

    // sys.mu = meV2au(20.0);

    // printer.printBandStructureSlice("data/BS.dat", Eigen::VectorXd::LinSpaced(1001, -0.5, 0.5), 0, 0.0);
    // printer.printBandStructure_discrete_ky("data/BS.dat", Eigen::VectorXd::LinSpaced(201, -0.3, 0.3), 15);
    // printer.printBandStructure_sparse_discrete_ky("data/BS.dat", Eigen::VectorXd::LinSpaced(201, -0.3, 0.3), 500, 80);


    // printer.printAbsDelta("data/abs_delta.dat", Eigen::VectorXd::LinSpaced(101, -0.5, 0.5), Eigen::VectorXd::LinSpaced(101, -0.5, 0.5));
    // printer.printAbelianBerryCurvature("data/BC.dat", Eigen::VectorXd::LinSpaced(101, -0.5, 0.5), Eigen::VectorXd::LinSpaced(101, -0.5, 0.5));
    // printer.printBerryCurvatureFromWilsonLoop("data/BC_WL.dat", Eigen::VectorXd::LinSpaced(101, -0.5, 0.5), Eigen::VectorXd::LinSpaced(101, -0.5, 0.5));
    // printer.printWilsonLoopSpectrum("data/WL.dat", 201, 5001, 0.3);
    // printer.printBandStructure_discrete_ky("data/BS.dat", Eigen::VectorXd::LinSpaced(201, -0.5, 0.5), 0.0, 251, true, 30, 0.0);

    // printer.printProbDen_sparse_discrete("data/prob_den.dat", 250, 9, meV2au(0.0));

    // discrete ky band structure
    {

        // sys.mu = meV2au(0.0);
        // sys.Bx = T2au(0.3);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(0.0);

        // std::size_t n_kx = 201;
        // double k_max = 0.3;
        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_kx, -k_max, k_max);

        // std::size_t n_ky = 251;

        // // std::cout << (sys.H_discrete_ky_invFourier(0.0, n_ky) - sys.H_discrete_ky(0.0, n_ky)).cwiseAbs().sum() / meV2au(1) << std::endl;

        // // printer.printBandStructure_discrete_ky("data/BS.dat", kx_vec, n_ky);
        // printer.printBandStructure_sparse_discrete_ky("data/BS.dat", kx_vec, n_ky, 80);
    }

    // discrete hamiltonian energy vs parameters
    {
        // System2D sys(HBdG, p);
        // sys.setHamiltonians_discrete(HBdG_discrete_onsite, HBdG_discrete_hopping_xp, HBdG_discrete_hopping_xm, HBdG_discrete_hopping_yp, HBdG_discrete_hopping_ym);

        // sys._p.mu = meV2au(4.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.0);

        // std::size_t n_kx = 250;
        // std::size_t n_ky = 9;
        // std::size_t n = 201;

        // // // kręcenie polem
        // //  double start = 0.0;
        // //  double end = 2.0 * M_PI;
        // //  Eigen::VectorXd vec = Eigen::VectorXd::LinSpaced(n, start, end);

        // // std::ofstream output_file("data/energy.dat");

        // // for (auto i = 0; i < n; ++i)
        // // {
        // //     sys._p.Bx = T2au(0.3) * std::cos(vec(i));
        // //     sys._p.By = T2au(0.3) * std::sin(vec(i));
        // //     // sys._p.mu = meV2au(vec(i));
        // //     auto [evals, evecs] = sys.H_discrete_eigendecomposiiton(n_kx, n_ky);
        // //     output_file << vec(i) << " " << evals.transpose() / meV2au(1) << std::endl;
        // // }

        // // samo pole
        // double start = 0.0;
        // double end = 0.5;
        // Eigen::VectorXd vec = Eigen::VectorXd::LinSpaced(n, start, end);

        // std::ofstream output_file("data/energy.dat");

        // for (auto i = 0; i < n; ++i)
        // {
        //     sys._p.Bz = T2au(vec(i));
        //     // sys._p.mu = meV2au(vec(i));
        //     auto [evals, evecs] = sys.H_discrete_eigendecomposiiton(n_kx, n_ky);
        //     // std::sort(evals.data(), evals.data() + evals.size());
        //     output_file << vec(i) << " " << evals.transpose() / meV2au(1) << std::endl;
        // }
    }

    // discrete system prob den
    {
        // System2D sys(HBdG, p);
        // sys.setHamiltonians_discrete(HBdG_discrete_onsite, HBdG_discrete_hopping_xp, HBdG_discrete_hopping_xm, HBdG_discrete_hopping_yp, HBdG_discrete_hopping_ym);

        // sys._p.mu = meV2au(4.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.3);
        // sys._p.Bz = T2au(1e-6);

        // std::size_t n_kx = 250;
        // std::size_t n_ky = 9;

        // double energy = meV2au(0.0);

        // sys.printProbDen_discrete("data/prob_den.dat", n_kx, n_ky, energy);
    }

    // delta vs momentum
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0); // doesnt matter for abs delta
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.2);

        // std::size_t n_k = 201;
        // double k_max = 0.3;

        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -k_max, k_max);
        // Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -k_max, k_max);

        // // sys.printAbsDelta("data/abs_delta.dat", kx_vec, ky_vec);
        // sys.printDeltaFromUnitaryTransformation("data/delta.dat", "data/DT.dat", kx_vec, ky_vec);
    }

    // abs delta at given k vs Bz
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0); // doesnt matter
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.0);

        // Point2D k{0.0, 0.0};
        // // Point2D k{M_PI, 0.0};
        // // Point2D k{0.0, M_PI};
        // // Point2D k{M_PI, M_PI};

        // std::size_t n_B = 1001;
        // double B_max = 5.0;
        // Eigen::VectorXd B_vec = Eigen::VectorXd::LinSpaced(n_B, -B_max, B_max);

        // std::ofstream output_file("data/abs_delta_B.dat");

        // output_file << "# k = " << k.transpose() << std::endl;

        // for (auto B : B_vec)
        // {
        //     sys._p.By = T2au(B);
        //     output_file << B << " " << sys.calcAbsDelta(k).transpose() / meV2au(1) << std::endl;
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

        // sys._p.mu = meV2au(0.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.3);

        // std::size_t n_dense = 100;  // number of k-points in each direction - even,
        // std::size_t n_sparse = 20; // number of k-points in each direction - even,

        // std::ofstream mu_Cs_file("data/CN.dat");

        // Eigen::VectorXd vec = Eigen::VectorXd::LinSpaced(51, -1.0, 1.0);

        // for (auto i = 0; i < vec.size(); ++i)
        // {
        //     // sys._p.Bz = T2au(vec(i));
        //     sys._p.mu = meV2au(vec(i));
        //     // mu_Cs_file << vec(i) << "\t" << sys.calcChernNumberUsingMatrixBerryCurvatureTrace(n_dense, n_sparse, 0.5) << std::endl;
        //     mu_Cs_file << vec(i) << "\t" << sys.calcCNUsingWilsonLoop(n_dense, n_sparse, 0.5) << std::endl;
        //     std::cout << i + 1 << " out of " << vec.size() << std::endl;
        // }
    }

    // Chern numbers vs parameters
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.0);

        // std::size_t n_dense = 16; // number of k-points in each direction - even,
        // std::size_t n_sparse = 6; // number of k-points in each direction - even,

        // Eigen::VectorXd p1_vec = Eigen::VectorXd::LinSpaced(201, -1, 5);
        // Eigen::VectorXd p2_vec = Eigen::VectorXd::LinSpaced(201, -1, 1);

        // std::ofstream params_CN_file("data/params_CN.dat");

        // std::size_t n_ky = 21;
        // sys.setHamiltonians_discrete_ky(HBdG_discrete_ky_onsite, HBdG_discrete_ky_hopping_p, HBdG_discrete_ky_hopping_m);

        // for (auto i = 0; i < p1_vec.size(); ++i)
        // {
        //     for (auto j = 0; j < p2_vec.size(); ++j)
        //     {
        //         sys._p.mu = meV2au(p1_vec(i));
        //         sys._p.Bx = T2au(p2_vec(j));
        //         // sys._p.By = T2au(p1_vec(i));;
        //         // sys._p.Bz = T2au(p2_vec(j));

        //         params_CN_file << p1_vec(i) << "\t"
        //                        << p2_vec(j) << "\t"
        //                        << sys.calcChernNumberUsingWilsonLoop_discrete_ky(n_dense, n_sparse, 0.7, n_ky) << "\n";
        //         //    << sys.calcCNUsingWilsonLoop(n_dense, n_sparse, 0.3) << std::endl;
        //         std::cout << i * p1_vec.size() + j + 1 << " out of " << p1_vec.size() * p2_vec.size() << std::endl;
        //     }
        // }
    }

    // BC
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0);
        // sys._p.Bx = T2au(0.0);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.2);

        // std::size_t n_k = 301;
        // double k_max = 0.5;

        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -k_max, k_max);
        // Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -k_max, k_max);

        // sys.printBerryCurvature("data/BC.dat", kx_vec, ky_vec);
        // sys.printBerryCurvatureFromWilsonLoop("data/BC_WL.dat", kx_vec, ky_vec);
    }

    // gap vs momentum
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0);
        // sys._p.Bx = T2au(0.3);
        // sys._p.By = T2au(0.0);
        // sys._p.Bz = T2au(0.0);

        // std::size_t n_k = 501;
        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -0.3, 0.3);
        // Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -0.3, 0.3);

        // sys.printDeltaFromUnitaryTransformation("data/delta.dat", "data/DT.dat", kx_vec, kx_vec);
    }

    // FS contours and gap along them
    {
        // sys.mu = meV2au(0.0);
        // sys.Bx = T2au(0.3);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(0.0);

        // auto FS = calc.FSContours();

        // for (auto ic = 0; ic < FS.size(); ++ic)
        // {
        //     printer.printAbsDeltaAlongContour("data/gap_FS" + std::to_string(ic) + ".dat", FS[ic]);
        // }

        // printer.printAbsDelta("data/abs_delta.dat", Eigen::VectorXd::LinSpaced(201, -0.3, 0.3), Eigen::VectorXd::LinSpaced(201, -0.3, 0.3));
    }

    // wilson loop spectrum
    {
        // System2D sys(HBdG, p);

        // sys._p.mu = meV2au(0.0);
        // sys._p.Bx = T2au(0.2);
        // sys._p.By = T2au(0.2);
        // sys._p.Bz = T2au(0.001);

        // sys.printWilsonLoopSpectrum("data/WL.dat", 501, 1001, 0.6);
    }

    return 0;
}