#include <iostream>
#include <iomanip>
#include <fstream>

#include "ToyModel.h"
#include "LAOSTO.h"

#include "System2DCalculations.h"
#include "System2DCalculationsPrinter.h"

int main()
{
    // ToyModel sys;
    LAOSTO sys;

    System2DCalculations calc(sys);
    System2DCalculationsPrinter printer(calc);

    double LAOSTO_bottom_band = -2.83022409703;
    double LAOSTO_top_band = 3.33333333333;

    // sys.mu = meV2au(LAOSTO_top_band);

    sys.mu = meV2au(4.0);
    sys.Bx = T2au(1.0);
    sys.By = T2au(1.0);
    sys.Bz = T2au(1.0);

    // std::cout << sys.HBdG_discrete(1,2).real()/meV2au(1) << std::endl; // Initialize the Hamiltonian

    
    // printer.printBandStructureSlice("data/BS_slice.dat", Eigen::VectorXd::LinSpaced(1001, -0.5, 0.5), 0);
    // printer.printBandStructure_discrete_ky("data/BS.dat", Eigen::VectorXd::LinSpaced(201, -0.5, 0.0), 80);
    // printer.printBandStructure_sparse_discrete_ky("data/BS.dat", Eigen::VectorXd::LinSpaced(501, -0.5, 0.0), 30, 2);

    // printer.printProbDen_sparse_discrete("data/prob_den.dat", 200, 9, 0.0);

    // for (size_t i = 1; i < 15; i++)
    // {
    //     /* code */
    //     // std::cout << i << " " << calc.ChernNumberUsingWilsonLoop_discrete_ky(50, 20, 0.4, i) << std::endl;
    // }

    // 2D k-space chern vs params
    {
        // sys.mu = meV2au(0.0);
        // sys.Bx = T2au(0.0);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(0.0);

        // // //TOYMODEL
        // // std::size_t n_par1 = 101;
        // // std::size_t n_par2 = 101;

        // // std::size_t n_dense = 20;
        // // std::size_t n_sparse = 10;

        // // Eigen::VectorXd vec_par1 = Eigen::VectorXd::LinSpaced(n_par1, -1.0, 1.0);
        // // Eigen::VectorXd vec_par2 = Eigen::VectorXd::LinSpaced(n_par2, -1.0, 1.0);

        // //LAOSTO
        // sys.mu = meV2au(LAOSTO_bottom_band);
        // std::size_t n_par1 = 51;
        // std::size_t n_par2 = 51;

        // std::size_t n_dense = 10;
        // std::size_t n_sparse = 5;

        // Eigen::VectorXd vec_par1 = Eigen::VectorXd::LinSpaced(n_par1, -4, 4);
        // Eigen::VectorXd vec_par2 = Eigen::VectorXd::LinSpaced(n_par2, -1, 1);

        // std::ofstream output_file("data/CN.dat");

        // for (auto i = 0; i < vec_par1.size(); ++i)
        // {
        //     for (auto j = 0; j < vec_par2.size(); ++j)
        //     {
        //         // sys.mu = meV2au(vec_par1(i));
        //         sys.Bx = T2au(vec_par1(i));
        //         // sys.By = T2au(vec_par1(i));
        //         // sys.Bz = T2au(vec_par1(i));

        //         // sys.Bx = T2au(vec_par2(j));
        //         // sys.By = T2au(vec_par2(j));
        //         sys.Bz = T2au(vec_par2(j));

        //         output_file << vec_par1(i) << " " << vec_par2(j) << " " << calc.ChernNumberUsingWilsonLoop(n_dense, n_sparse, 0.5) << std::endl;
        //         // output_file << vec_par1(i) << " " << vec_par2(j) << " " << calc.ChernNumberUsingBerryCurvatureFromWilsonLoop(n_dense, n_sparse, 0.5) << std::endl;
        //     }
        // }
    }

    // 1D k-space chern vs params
    {
        // sys.mu = meV2au(0.0);
        // sys.Bx = T2au(0.0);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(0.0);

        // // // TOYMODEL
        // // std::size_t n_ky_min = 1;
        // // std::size_t n_ky_max = 251;

        // // std::size_t n_par = 401;
        // // Eigen::VectorXd vec_par = Eigen::VectorXd::LinSpaced(n_par, -2.0, 2.0);

        // // std::size_t n_dense = 20;
        // // std::size_t n_sparse = 10;

        // // //LAOSTO
        // sys.mu = meV2au(LAOSTO_top_band);
        // std::size_t n_ky_min = 1;
        // std::size_t n_ky_max = 30;

        // std::size_t n_par = 51;
        // Eigen::VectorXd vec_par = Eigen::VectorXd::LinSpaced(n_par, 0.0, 4.0);

        // std::size_t n_dense = 10;
        // std::size_t n_sparse = 5;

        // std::ofstream output_file("data/CN1D.dat");

        // for (auto n_ky = n_ky_min; n_ky <= n_ky_max; ++n_ky)
        // {
        //     for (auto j = 0; j < vec_par.size(); ++j)
        //     {
        //         sys.Bx = T2au(vec_par(j));
        //         // sys.By = T2au(vec_par(j));
        //         // sys.Bz = T2au(vec_par(j));

        //         output_file << n_ky << " " << vec_par(j) << " " << calc.ChernNumberUsingWilsonLoop_discrete_ky(n_dense, n_sparse, 0.5, n_ky) << std::endl;
        //         // output_file << vec_par1(i) << " " << vec_par2(j) << " " << calc.ChernNumberUsingBerryCurvatureFromWilsonLoop(n_dense, n_sparse, 0.5) << std::endl;
        //     }
        // }
    }

    // discrete hamiltonian energy vs parameters
    {

        // sys.mu = meV2au(4.0);
        // sys.Bx = T2au(0.0);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(0.0);

        // std::size_t n_kx = 200;
        // std::size_t n_ky = 9;
        // std::size_t n = 101;

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
        // double end = 1.0;
        // Eigen::VectorXd vec = Eigen::VectorXd::LinSpaced(n, start, end);

        // std::ofstream output_file("data/energy.dat");

        // for (auto i = 0; i < n; ++i)
        // {
        //     sys.Bz = T2au(vec(i));
        //     // sys.mu = meV2au(vec(i));
        //     auto evals = calc.eigenvals_sparse_discrete(n_kx, n_ky, 30);
        //     std::sort(evals.data(), evals.data() + evals.size());
        //     output_file << vec(i) << " " << evals.transpose() / meV2au(1) << std::endl;
        // }
    }

    // discrete system prob den
    {
        // sys.mu = meV2au(4.0);
        // sys.Bx = T2au(1.3);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(0.0);

        // std::size_t n_kx = 250;
        // std::size_t n_ky = 9;

        // double energy = meV2au(0.0);

        // printer.printProbDen_sparse_discrete("data/prob_den.dat", n_kx, n_ky, energy);
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
        // sys.mu = meV2au(0.0);
        // sys.Bx = T2au(0.0);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(0.0);

        // std::size_t n_dense = 16; // number of k-points in each direction - even,
        // std::size_t n_sparse = 6; // number of k-points in each direction - even,

        // Eigen::VectorXd p1_vec = Eigen::VectorXd::LinSpaced(13, 0, 6);
        // Eigen::VectorXd p2_vec = Eigen::VectorXd::LinSpaced(17, 2, 10);

        // std::ofstream params_CN_file("data/params_CN.dat");

        // std::size_t n_ky = 21;

        // for (auto i = 0; i < p1_vec.size(); ++i)
        // {
        //     for (auto j = 0; j < p2_vec.size(); ++j)
        //     {
        //         sys.mu = meV2au(p1_vec(i));
        //         // sys.Bx = T2au(p2_vec(j));
        //         // sys.By = T2au(p1_vec(i));;
        //         sys.Bz = T2au(p2_vec(j));

        //         params_CN_file << p1_vec(i) << "\t"
        //                        << p2_vec(j) << "\t"
        //                        << calc.ChernNumberUsingWilsonLoop_discrete_ky(n_dense, n_sparse, 0.5, n_ky) << std::endl;
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
        // sys.mu = meV2au(0.0);
        // sys.Bx = T2au(2.0);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(0.0);

        // printer.printWilsonLoopSpectrum("data/WL.dat", 101, 5001, 0.3);
    }

    return 0;
}