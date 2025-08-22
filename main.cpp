#include <iostream>
#include <iomanip>
#include <fstream>

#include "ToyModel.h"
#include "LAOSTO.h"

#include "System2DCalculations.h"
#include "System2DCalculationsPrinter.h"

int main(int argc, char *argv[])
{
    // ToyModel sys;
    LAOSTO sys;

    System2DCalculations calc(sys);
    System2DCalculationsPrinter printer(calc);

    double LAOSTO_low_band = -47.5031092363;
    double LAOSTO_bottom_band = -2.83022409703;
    double LAOSTO_top_band = 3.33333333333;

    // std::cout << std::setprecision(12) << calc.eigenvals_discrete_ky_normal(0.0, 12)/meV2au(1) << std::endl;

    // sys.mu = meV2au(19.5773907141);
    // sys.Bx = T2au(0.184);

    // printer.printProbDen_sparse_discrete("data/prob_den.dat", 5000, 12, 0.0);

    // // //LAOSTO MZM nky=12, nkx=5000
    // sys.mu = meV2au(-1.55022);
    // sys.Bx = T2au(1.0);
    // sys.By = T2au(0.0);
    // sys.Bz = T2au(0.0);

    // sys.mu = meV2au(11.30);
    // sys.Bx = T2au(0.5);
    // sys.By = T2au(0.0);
    // sys.Bz = T2au(0.0);
    // std::cout << sys.HBdG_discrete(1,2).real()/meV2au(1) << std::endl; // Initialize the Hamiltonian

    // printer.printBandStructureSlice_normal("data/BS_slice.dat", Eigen::VectorXd::LinSpaced(1001, -0.5, 0.5), 0);
    // printer.printBandStructureSlice("data/BS_slice.dat", Eigen::VectorXd::LinSpaced(1001, -0.5, 0.5), 0);
    // printer.printBandStructure_discrete_ky("data/BS.dat", Eigen::VectorXd::LinSpaced(201, -0.5, 0.0), 80);
    // printer.printBandStructure_sparse_discrete_ky("data/BS.dat", Eigen::VectorXd::LinSpaced(201, -0.2, 0.2), 200, 30);

    // printer.printBerryCurvatureFromWilsonLoop("data/BC.dat", Eigen::VectorXd::LinSpaced(201, -0.5, 0.5), Eigen::VectorXd::LinSpaced(201, -0.5, 0.5));

    // printer.printProbDen_sparse_discrete("data/prob_den.dat", 5000, 12, 0.0);

    // for (size_t i = 1; i < 15; i++)
    // {
    //     /* code */
    //     // std::cout << i << " " << calc.ChernNumberUsingWilsonLoop_discrete_ky(50, 20, 0.4, i) << std::endl;
    // }

    // Z2 Pfaffian testing - 2D case
    {
        Eigen::MatrixXd I_n_bands = Eigen::MatrixXd::Identity(sys.n_bands, sys.n_bands);
        Eigen::MatrixXcd U(2 * sys.n_bands, 2 * sys.n_bands);

        U << I_n_bands, I_n_bands,
            -1i * I_n_bands, 1i * I_n_bands;

        for (auto mu : Eigen::VectorXd::LinSpaced(101, 3, 4))
        {

            for (auto Bz : Eigen::VectorXd::LinSpaced(101, -2, 2))
            {
                sys.mu = meV2au(mu);
                sys.Bx = T2au(0.0);
                sys.By = T2au(0.0);
                sys.Bz = T2au(Bz);

                Eigen::MatrixXcd H1 = sys.HBdG(0.0, 0.0);
                Eigen::MatrixXcd H2 = sys.HBdG(M_PI, 0.0);
                Eigen::MatrixXcd H3 = sys.HBdG(0.0, M_PI);
                Eigen::MatrixXcd H4 = sys.HBdG(M_PI, M_PI);

                Eigen::MatrixXd A1 = (U * H1 * U.adjoint()).imag();
                Eigen::MatrixXd A2 = (U * H2 * U.adjoint()).imag();
                Eigen::MatrixXd A3 = (U * H3 * U.adjoint()).imag();
                Eigen::MatrixXd A4 = (U * H4 * U.adjoint()).imag();

                double pf1 = pfaffian(A1);
                double pf2 = pfaffian(A2);
                double pf3 = pfaffian(A3);
                double pf4 = pfaffian(A4);

                std::cout << mu << " " << Bz << " " << 0.5 - 0.5 * (pf1 * pf2 * pf3 * pf4 > 0.0 ? 1.0 : -1.0) << std::endl;
            }
        }
    }

    // Z2 Pfaffian testing - 1D case
    {
        // Eigen::MatrixXd I_n_bands = Eigen::MatrixXd::Identity(sys.n_bands, sys.n_bands);
        // Eigen::MatrixXcd U_single(2 * sys.n_bands, 2 * sys.n_bands);

        // U_single << I_n_bands, I_n_bands,
        //     -1i * I_n_bands, 1i * I_n_bands;

        // std::size_t n_ky = 12;
        // Eigen::MatrixXcd U = Eigen::kroneckerProduct(Eigen::MatrixXcd::Identity(n_ky, n_ky), U_single).eval();

        // for (auto mu : Eigen::VectorXd::LinSpaced(51, 4.4, 5.0))
        // {

        //     for (auto Bx : Eigen::VectorXd::LinSpaced(51, -2, 2))
        //     {
        //         sys.mu = meV2au(mu);
        //         sys.Bx = T2au(Bx);
        //         sys.By = T2au(0.0);
        //         sys.Bz = T2au(0.0);

        //         Eigen::MatrixXcd H1 = sys.HBdG_discrete_ky(0.0, n_ky);
        //         Eigen::MatrixXcd H2 = sys.HBdG_discrete_ky(M_PI, n_ky);

        //         Eigen::MatrixXd A1 = (U * H1 * U.adjoint()).imag();
        //         Eigen::MatrixXd A2 = (U * H2 * U.adjoint()).imag();

        //         double pf1 = pfaffian(A1);
        //         double pf2 = pfaffian(A2);

        //         // double det1 = A1.determinant();
        //         // double det2 = A2.determinant();

        //         // std::cout << "Pf1 = " << pf1 << std::endl;
        //         // std::cout << "det1 = " << det1 << std::endl;
        //         // std::cout << "|det1-pf1^2| = " << std::abs(pf1 * pf1 - det1) << std::endl;

        //         // std::cout << "Pf2 = " << pf2 << std::endl;
        //         // std::cout << "det2 = " << det2 << std::endl;
        //         // std::cout << "|det2-pf2^2| = " << std::abs(pf2 * pf2 - det2) << std::endl;

        //         std::cout << mu << " " << Bx << " " << 0.5 - 0.5 * (pf1 * pf2 > 0.0 ? 1.0 : -1.0) << "\n";
        //     }
        // }
    }

    // 2D k-space chern vs params
    {
        // sys.mu = meV2au(0.0);
        // sys.Bx = T2au(0.0);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(0.0);

        // // //TOYMODEL
        // // std::size_t n_par1 = 51;
        // // std::size_t n_par2 = 51;

        // // std::size_t n_dense = 10;
        // // std::size_t n_sparse = 5;

        // // Eigen::VectorXd vec_par1 = Eigen::VectorXd::LinSpaced(n_par1, -1.0, 1.0);
        // // Eigen::VectorXd vec_par2 = Eigen::VectorXd::LinSpaced(n_par2, -1.0, 1.0);

        // //LAOSTO
        // // sys.mu = meV2au(LAOSTO_top_band);

        // std::size_t n_par1 = 201;
        // std::size_t n_par2 = 201;

        // std::size_t n_dense = 16;
        // std::size_t n_sparse = 8;

        // Eigen::VectorXd vec_par1 = Eigen::VectorXd::LinSpaced(n_par1, -50, -40);
        // Eigen::VectorXd vec_par2 = Eigen::VectorXd::LinSpaced(n_par2, -5, 5);

        // std::ofstream output_file("data/CN.dat");

        // for (auto i = 0; i < vec_par1.size(); ++i)
        // {
        //     for (auto j = 0; j < vec_par2.size(); ++j)
        //     {
        //         sys.mu = meV2au(vec_par1(i));
        //         // sys.Bx = T2au(vec_par1(i));
        //         // sys.By = T2au(vec_par1(i));
        //         // sys.Bz = T2au(vec_par1(i));

        //         // sys.Bx = T2au(vec_par2(j));
        //         sys.By = T2au(vec_par2(j));
        //         // sys.Bz = T2au(vec_par2(j));

        //         output_file << vec_par1(i) << " " << vec_par2(j) << " " << calc.ChernNumberUsingWilsonLoop(n_dense, n_sparse, 0.5) << std::endl;
        //         // output_file << vec_par1(i) << " " << vec_par2(j) << " " << calc.ChernNumberUsingBerryCurvatureFromWilsonLoop(n_dense, n_sparse, 0.5) << std::endl;
        //     }
        // }
    }

    // 2D k-space chern vs single parameter
    {
        // sys.mu = meV2au(0.0);
        // sys.Bx = T2au(0.0);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(5.0);

        // std::size_t n_par = 501;

        // std::size_t n_dense = 16;
        // std::size_t n_sparse = 8;

        // Eigen::VectorXd vec_par = Eigen::VectorXd::LinSpaced(n_par, -55, -45);

        // std::ofstream output_file("data/CN_single_param.dat");

        // for (auto i = 0; i < vec_par.size(); ++i)
        // {
        //         sys.mu = meV2au(vec_par(i));
        //         // sys.Bx = T2au(vec_par(i));
        //         // sys.By = T2au(vec_par(i));
        //         // sys.Bz = T2au(vec_par(i));

        //         output_file << vec_par(i) << " " << calc.ChernNumberUsingWilsonLoop(n_dense, n_sparse, 0.5) << std::endl;
        // }
    }

    // 1D k-space chern vs params
    {
        // sys.mu = meV2au(0.0);
        // sys.Bx = T2au(0.0);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(0.0);

        // // //TOYMODEL
        // // std::size_t n_par1 = 51;
        // // std::size_t n_par2 = 51;

        // // std::size_t n_dense = 10;
        // // std::size_t n_sparse = 5;

        // // Eigen::VectorXd vec_par1 = Eigen::VectorXd::LinSpaced(n_par1, -1.0, 1.0);
        // // Eigen::VectorXd vec_par2 = Eigen::VectorXd::LinSpaced(n_par2, -1.0, 1.0);

        // // LAOSTO
        // std::size_t n_ky = 12;

        // std::size_t n_par1 = 401;
        // std::size_t n_par2 = 401;

        // std::size_t n_dense = 16;
        // std::size_t n_sparse = 8;

        // Eigen::VectorXd vec_par1 = Eigen::VectorXd::LinSpaced(n_par1, -2, 2);
        // Eigen::VectorXd vec_par2 = Eigen::VectorXd::LinSpaced(n_par2, -2, 2);

        // std::ofstream output_file("data/CN1D.dat");

        // for (auto i = 0; i < vec_par1.size(); ++i)
        // {
        //     for (auto j = 0; j < vec_par2.size(); ++j)
        //     {
        //         sys.mu = meV2au(vec_par1(i));
        //         // sys.Bx = T2au(vec_par1(i));
        //         // sys.By = T2au(vec_par1(i));
        //         // sys.Bz = T2au(vec_par1(i));

        //         sys.Bx = T2au(vec_par2(j));
        //         // sys.By = T2au(vec_par2(j));
        //         // sys.Bz = T2au(vec_par2(j));

        //         // output_file << vec_par1(i) << " " << vec_par2(j) << " " << calc.ChernNumberUsingWilsonLoop(n_dense, n_sparse, 0.5) << std::endl;
        // output_file << vec_par1(i) << " " << vec_par2(j) << " " << calc.ChernNumberUsingWilsonLoop_discrete_ky(n_dense, n_sparse, 0.5, n_ky) << std::endl;
        //         // output_file << vec_par1(i) << " " << vec_par2(j) << " " << calc.ChernNumberUsingBerryCurvatureFromWilsonLoop(n_dense, n_sparse, 0.5) << std::endl;
        //     }
        // }
    }

    // 1D k-space chern vs nky and param
    {
        //     sys.mu = meV2au(0.0);
        //     sys.Bx = T2au(1.0);
        //     sys.By = T2au(0.0);
        //     sys.Bz = T2au(0.0);

        //     // // TOYMODEL
        //     // std::size_t n_ky_min = 1;
        //     // std::size_t n_ky_max = 40;

        //     // std::size_t n_par = 51;
        //     // Eigen::VectorXd vec_par = Eigen::VectorXd::LinSpaced(n_par, -2.0, 2.0);

        //     // std::size_t n_dense = 10;
        //     // std::size_t n_sparse = 5;

        //     // //LAOSTO
        // //    12 -1.55022 -0.5
        //     std::size_t n_ky_min = 12;
        //     std::size_t n_ky_max = 12;

        //     std::size_t n_par = 101;
        //     Eigen::VectorXd vec_par = Eigen::VectorXd::LinSpaced(n_par, LAOSTO_bottom_band+1, LAOSTO_bottom_band+2);

        //     std::size_t n_dense = 16;
        //     std::size_t n_sparse = 5;

        //     std::ofstream output_file("data/CN1D.dat");

        //     for (auto n_ky = n_ky_min; n_ky <= n_ky_max; ++n_ky)
        //     {
        //         for (auto j = 0; j < vec_par.size(); ++j)
        //         {
        //             sys.mu = meV2au(vec_par(j));
        //             // sys.Bx = T2au(vec_par(j));
        //             // sys.By = T2au(vec_par(j));
        //             // sys.Bz = T2au(vec_par(j));

        //             output_file << n_ky << " " << vec_par(j) << " " << calc.ChernNumberUsingWilsonLoop_discrete_ky(n_dense, n_sparse, 0.5, n_ky) << std::endl;
        //             // output_file << vec_par1(i) << " " << vec_par2(j) << " " << calc.ChernNumberUsingBerryCurvatureFromWilsonLoop(n_dense, n_sparse, 0.5) << std::endl;
        //         }
        //     }
    }

    // 1D k-space chern vs single param
    {
        // sys.mu = meV2au(0);
        // sys.Bx = T2au(1.0);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(0.0);

        // // LAOSTO
        // std::size_t n_ky = 12;

        // std::size_t n_par1 = 40*100+1;

        // std::size_t n_dense = 16;
        // std::size_t n_sparse = 8;

        // Eigen::VectorXd vec_par1 = Eigen::VectorXd::LinSpaced(n_par1, 0, 40);

        // std::ofstream output_file("data/CN1D_single_param.dat");

        // for (auto i = 0; i < vec_par1.size(); ++i)
        // {
        //     sys.mu = meV2au(vec_par1(i));
        //     // sys.Bx = T2au(vec_par1(i));
        //     // sys.By = T2au(vec_par1(i));
        //     // sys.Bz = T2au(vec_par1(i));

        //     output_file << vec_par1(i) << " " << calc.ChernNumberUsingWilsonLoop_discrete_ky(n_dense, n_sparse, 0.7, n_ky) << std::endl;
        // }
    }

    // discrete hamiltonian energy vs parameters
    {

        // sys.mu = meV2au(9.53);
        // sys.Bx = T2au(0.0);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(0.0);

        // std::size_t n_kx = 500;
        // std::size_t n_ky = 12;
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

        // double start = 0;
        // double end = 2;
        // Eigen::VectorXd vec = Eigen::VectorXd::LinSpaced(n, start, end);

        // std::ofstream output_file("data/energy.dat");

        // for (auto i = 0; i < n; ++i)
        // {
        //     // sys.mu = meV2au(vec(i));
        //     sys.Bx = T2au(vec(i));
        //     // sys.Bz = T2au(vec(i));

        //     auto evals = calc.eigenvals_sparse_discrete(n_kx, n_ky, 80);
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
        // sys.mu = meV2au(0);
        // sys.Bx = T2au(0.0);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(-1.0);

        // printer.printWilsonLoopSpectrum("data/WL.dat", 101, 2001, 0.4);
    }

    return 0;
}