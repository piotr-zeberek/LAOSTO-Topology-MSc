#include <iostream>
#include <iomanip>
#include <fstream>

#include "ToyModel.h"
#include "LAOSTO.h"
#include "LAOSTO_exts.h"

#include "System2DCalculations.h"
#include "System2DCalculationsPrinter.h"

int main(int argc, char *argv[])
{
    // ToyModel sys;
    // LAOSTO sys;
    LAOSTO_exts sys;

    System2DCalculations calc(sys);
    System2DCalculationsPrinter printer(calc);

    double LAOSTO_low_band = -47.5031092363;
    double LAOSTO_bottom_band = -2.83022409703;
    double LAOSTO_top_band = 3.33333333333;

    sys.mu = meV2au(LAOSTO_bottom_band);
    sys.Bx = T2au(0.0);
    sys.By = T2au(0.0);
    sys.Bz = T2au(1);

    printer.printBandStructureSlice_orbital_type("data/BS_slice_orbital.dat", Eigen::VectorXd::LinSpaced(2001, -0.1, 0.1), 0);


    // 2D k-space Z2 from Pfaffian vs 2 params - CLI arguments
    // Z2 init_mu init_Bx init_By init_Bz kind_par1 n_par1 min_par1 max_par1 kind_par2 n_par2 min_par2 max_par2 output_filename
    // Parameters:
    // 0: program name
    // 1: init_mu
    // 2: init_Bx
    // 3: init_By
    // 4: init_Bz
    // 5: kind_par1 (0 - mu, 1 - Bx, 2 - By, 3 - Bz)
    // 6: n_par1
    // 7: min_par1
    // 8: max_par1
    // 9: kind_par2 (0 - mu, 1 - Bx, 2 - By, 3 - Bz)
    // 10: n_par2
    // 11: min_par2
    // 12: max_par2
    // 13: output_filename
    {
        // Eigen::MatrixXd I_n_bands = Eigen::MatrixXd::Identity(sys.n_bands, sys.n_bands);
        // Eigen::MatrixXcd U(2 * sys.n_bands, 2 * sys.n_bands);

        // U << I_n_bands, I_n_bands,
        //     -1i * I_n_bands, 1i * I_n_bands;

        // if (argc < 14)
        // {
        //     std::cerr << "Usage: " << argv[0]
        //               << " init_mu init_Bx init_By init_Bz kind_par1 n_par1 min_par1 max_par1 kind_par2 n_par2 min_par2 max_par2 output_filename"
        //               << std::endl;
        //     return 1;
        // }
        // sys.mu = meV2au(std::stod(argv[1]));
        // sys.Bx = T2au(std::stod(argv[2]));
        // sys.By = T2au(std::stod(argv[3]));
        // sys.Bz = T2au(std::stod(argv[4]));

        // std::size_t kind_par1 = std::stoul(argv[5]);
        // std::size_t n_par1 = std::stoul(argv[6]);
        // double min_par1 = std::stod(argv[7]);
        // double max_par1 = std::stod(argv[8]);

        // std::size_t kind_par2 = std::stoul(argv[9]);
        // std::size_t n_par2 = std::stoul(argv[10]);
        // double min_par2 = std::stod(argv[11]);
        // double max_par2 = std::stod(argv[12]);

        // std::string output_filename = argv[13];

        // Eigen::VectorXd vec_par1 = Eigen::VectorXd::LinSpaced(n_par1, min_par1, max_par1);
        // Eigen::VectorXd vec_par2 = Eigen::VectorXd::LinSpaced(n_par2, min_par2, max_par2);

        // std::ofstream output_file(output_filename);

        // for (auto i = 0; i < n_par1; ++i)
        // {
        //     for (auto j = 0; j < n_par2; ++j)
        //     {
        //         switch (kind_par1)
        //         {
        //         case 0: // mu
        //             sys.mu = meV2au(vec_par1(i));
        //             break;
        //         case 1: // Bx
        //             sys.Bx = T2au(vec_par1(i));
        //             break;
        //         case 2: // By
        //             sys.By = T2au(vec_par1(i));
        //             break;
        //         case 3: // Bz
        //             sys.Bz = T2au(vec_par1(i));
        //             break;
        //         default:
        //             std::cerr << "Invalid kind_par1: " << kind_par1 << std::endl;
        //             return 1;
        //         }

        //         switch (kind_par2)
        //         {
        //         case 0: // mu
        //             sys.mu = meV2au(vec_par2(j));
        //             break;
        //         case 1: // Bx
        //             sys.Bx = T2au(vec_par2(j));
        //             break;
        //         case 2: // By
        //             sys.By = T2au(vec_par2(j));
        //             break;
        //         case 3: // Bz
        //             sys.Bz = T2au(vec_par2(j));
        //             break;
        //         default:
        //             std::cerr << "Invalid kind_par2: " << kind_par2 << std::endl;
        //             return 1;
        //         }

        //         Eigen::MatrixXcd H1 = sys.HBdG(0.0, 0.0);
        //         Eigen::MatrixXcd H2 = sys.HBdG(0.0, M_PI);
        //         Eigen::MatrixXcd H3 = sys.HBdG(M_PI, M_PI);
        //         Eigen::MatrixXcd H4 = sys.HBdG(M_PI, 0.0);

        //         Eigen::MatrixXd A1 = (U * H1 * U.adjoint()).imag();
        //         Eigen::MatrixXd A2 = (U * H2 * U.adjoint()).imag();
        //         Eigen::MatrixXd A3 = (U * H3 * U.adjoint()).imag();
        //         Eigen::MatrixXd A4 = (U * H4 * U.adjoint()).imag();

        //         double pf1 = pfaffian(A1);
        //         double pf2 = pfaffian(A2);
        //         double pf3 = pfaffian(A3);
        //         double pf4 = pfaffian(A4);

        //         auto Z2 = [](double pfaffian_product) -> double
        //         { return (pfaffian_product > 0.0 ? 1.0 : -1.0); };

        //         output_file << vec_par1(i) << " "
        //                     << vec_par2(j) << " "
        //                     << Z2(pf1 * pf2 * pf3 * pf4) << " "
        //                     << Z2(pf1 * pf2) << " "
        //                     << Z2(pf2 * pf3) << " "
        //                     << Z2(pf3 * pf4) << " "
        //                     << Z2(pf4 * pf1) << "\n";
        //     }
        // }
    }

    // 1D k-space Z2 from Pfaffian vs 2 params - CLI arguments
    // Z2 init_mu init_Bx init_By init_Bz kind_par1 n_par1 min_par1 max_par1 kind_par2 n_par2 min_par2 max_par2 n_ky output_filename
    // Parameters:
    // 0: program name
    // 1: init_mu
    // 2: init_Bx
    // 3: init_By
    // 4: init_Bz
    // 5: kind_par1 (0 - mu, 1 - Bx, 2 - By, 3 - Bz)
    // 6: n_par1
    // 7: min_par1
    // 8: max_par1
    // 9: kind_par2 (0 - mu, 1 - Bx, 2 - By, 3 - Bz)
    // 10: n_par2
    // 11: min_par2
    // 12: max_par2
    // 13: n_ky
    // 14: output_filename
    {
        // if (argc < 15)
        // {
        //     std::cerr << "Usage: " << argv[0]
        //               << " init_mu init_Bx init_By init_Bz kind_par1 n_par1 min_par1 max_par1 kind_par2 n_par2 min_par2 max_par2 n_ky output_filename"
        //               << std::endl;
        //     return 1;
        // }
        // sys.mu = meV2au(std::stod(argv[1]));
        // sys.Bx = T2au(std::stod(argv[2]));
        // sys.By = T2au(std::stod(argv[3]));
        // sys.Bz = T2au(std::stod(argv[4]));

        // std::size_t kind_par1 = std::stoul(argv[5]);
        // std::size_t n_par1 = std::stoul(argv[6]);
        // double min_par1 = std::stod(argv[7]);
        // double max_par1 = std::stod(argv[8]);

        // std::size_t kind_par2 = std::stoul(argv[9]);
        // std::size_t n_par2 = std::stoul(argv[10]);
        // double min_par2 = std::stod(argv[11]);
        // double max_par2 = std::stod(argv[12]);

        // std::size_t n_ky = std::stoul(argv[13]);

        // std::string output_filename = argv[14];

        // Eigen::MatrixXd I_n_bands = Eigen::MatrixXd::Identity(sys.n_bands, sys.n_bands);
        // Eigen::MatrixXcd U_single(2 * sys.n_bands, 2 * sys.n_bands);

        // U_single << I_n_bands, I_n_bands,
        //     -1i * I_n_bands, 1i * I_n_bands;

        // Eigen::MatrixXcd U = Eigen::kroneckerProduct(Eigen::MatrixXcd::Identity(n_ky, n_ky), U_single).eval();

        // Eigen::VectorXd vec_par1 = Eigen::VectorXd::LinSpaced(n_par1, min_par1, max_par1);
        // Eigen::VectorXd vec_par2 = Eigen::VectorXd::LinSpaced(n_par2, min_par2, max_par2);

        // std::ofstream output_file(output_filename);

        // for (auto i = 0; i < n_par1; ++i)
        // {
        //     for (auto j = 0; j < n_par2; ++j)
        //     {
        //         switch (kind_par1)
        //         {
        //         case 0: // mu
        //             sys.mu = meV2au(vec_par1(i));
        //             break;
        //         case 1: // Bx
        //             sys.Bx = T2au(vec_par1(i));
        //             break;
        //         case 2: // By
        //             sys.By = T2au(vec_par1(i));
        //             break;
        //         case 3: // Bz
        //             sys.Bz = T2au(vec_par1(i));
        //             break;
        //         default:
        //             std::cerr << "Invalid kind_par1: " << kind_par1 << std::endl;
        //             return 1;
        //         }

        //         switch (kind_par2)
        //         {
        //         case 0: // mu
        //             sys.mu = meV2au(vec_par2(j));
        //             break;
        //         case 1: // Bx
        //             sys.Bx = T2au(vec_par2(j));
        //             break;
        //         case 2: // By
        //             sys.By = T2au(vec_par2(j));
        //             break;
        //         case 3: // Bz
        //             sys.Bz = T2au(vec_par2(j));
        //             break;
        //         default:
        //             std::cerr << "Invalid kind_par2: " << kind_par2 << std::endl;
        //             return 1;
        //         }

        //         Eigen::MatrixXcd H1 = sys.HBdG_discrete_ky(0.0, n_ky);
        //         Eigen::MatrixXcd H2 = sys.HBdG_discrete_ky(M_PI, n_ky);

        //         Eigen::MatrixXd A1 = (U * H1 * U.adjoint()).imag();
        //         Eigen::MatrixXd A2 = (U * H2 * U.adjoint()).imag();

        //         double pf1 = pfaffian(A1);
        //         double pf2 = pfaffian(A2);

        //         auto Z2 = [](double pfaffian_product) -> double
        //         { return (pfaffian_product > 0.0 ? 1.0 : -1.0); };

        //         output_file << vec_par1(i) << " "
        //                     << vec_par2(j) << " "
        //                     << Z2(pf1 * pf2) << "\n";
        //     }
        // }
    }

    // Z2 Pfaffian testing - 2D case
    {
        // Eigen::MatrixXd I_n_bands = Eigen::MatrixXd::Identity(sys.n_bands, sys.n_bands);
        // Eigen::MatrixXcd U(2 * sys.n_bands, 2 * sys.n_bands);

        // U << I_n_bands, I_n_bands,
        //     -1i * I_n_bands, 1i * I_n_bands;

        // for (auto mu : Eigen::VectorXd::LinSpaced(101, -1, 1))
        // {

        //     for (auto Bz : Eigen::VectorXd::LinSpaced(101, -1, 1))
        //     {
        //         sys.mu = meV2au(mu);
        //         sys.Bx = T2au(0.0);
        //         sys.By = T2au(0.0);
        //         sys.Bz = T2au(Bz);

        //         Eigen::MatrixXcd H1 = sys.HBdG(0.0, 0.0);
        //         Eigen::MatrixXcd H2 = sys.HBdG(M_PI, 0.0);
        //         Eigen::MatrixXcd H3 = sys.HBdG(0.0, M_PI);
        //         Eigen::MatrixXcd H4 = sys.HBdG(M_PI, M_PI);

        //         Eigen::MatrixXd A1 = (U * H1 * U.adjoint()).imag();
        //         Eigen::MatrixXd A2 = (U * H2 * U.adjoint()).imag();
        //         Eigen::MatrixXd A3 = (U * H3 * U.adjoint()).imag();
        //         Eigen::MatrixXd A4 = (U * H4 * U.adjoint()).imag();

        //         double pf1 = pfaffian(A1);
        //         double pf2 = pfaffian(A2);
        //         double pf3 = pfaffian(A3);
        //         double pf4 = pfaffian(A4);

        //         auto Z2 = [](double pfaffian_product) -> double
        //         { return (pfaffian_product > 0.0 ? 1.0 : -1.0); };

        //         std::cout << mu << " "
        //                   << Bz << " "
        //                   << Z2(pf1 * pf2 * pf3 * pf4) << " "
        //                   << Z2(pf1 * pf2) << " "
        //                   << Z2(pf2 * pf3) << " "
        //                   << Z2(pf3 * pf4) << " "
        //                   << Z2(pf4 * pf1) << "\n";
        //     }
        // }
    }

    // Z2 Pfaffian testing - 1D case
    {
        // Eigen::MatrixXd I_n_bands = Eigen::MatrixXd::Identity(sys.n_bands, sys.n_bands);
        // Eigen::MatrixXcd U_single(2 * sys.n_bands, 2 * sys.n_bands);

        // U_single << I_n_bands, I_n_bands,
        //     -1i * I_n_bands, 1i * I_n_bands;

        // std::size_t n_ky = 80;
        // Eigen::MatrixXcd U = Eigen::kroneckerProduct(Eigen::MatrixXcd::Identity(n_ky, n_ky), U_single).eval();

        // for (auto mu : Eigen::VectorXd::LinSpaced(1, 0, 0))
        // {

        //     for (auto B : Eigen::VectorXd::LinSpaced(401, -2, 2))
        //     {
        //         sys.mu = meV2au(mu);
        //         sys.Bx = T2au(B);
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

        //         std::cout << mu << " " << B << " " << (pf1 > 0.0 ? 1.0 : -1.0) * (pf2 > 0.0 ? 1.0 : -1.0) << "\n";
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

        // std::size_t n_par1 = 101;
        // std::size_t n_par2 = 101;

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
        //         // sys.By = T2au(vec_par2(j));
        //         sys.Bz = T2au(vec_par2(j));

        //         output_file << vec_par1(i) << " " << vec_par2(j) << " " << calc.ChernNumberUsingWilsonLoop(n_dense, n_sparse, 0.5) << std::endl;
        //         // output_file << vec_par1(i) << " " << vec_par2(j) << " " << calc.ChernNumberUsingBerryCurvatureFromWilsonLoop(n_dense, n_sparse, 0.5) << std::endl;
        //     }
        // }
    }

    // 2D k-space chern vs single parameter
    {
        // sys.mu = meV2au(LAOSTO_bottom_band);
        // sys.Bx = T2au(0.0);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(0.0);

        // std::size_t n_par = 501;

        // std::size_t n_dense = 16;
        // std::size_t n_sparse = 8;

        // Eigen::VectorXd vec_par = Eigen::VectorXd::LinSpaced(n_par, -5, 5);

        // std::ofstream output_file("data/CN_single_param.dat");

        // for (auto i = 0; i < vec_par.size(); ++i)
        // {
        //         // sys.mu = meV2au(vec_par(i));
        //         // sys.Bx = T2au(vec_par(i));
        //         // sys.By = T2au(vec_par(i));
        //         sys.Bz = T2au(vec_par(i));

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
        // sys.mu = meV2au(4.68);
        // sys.Bx = T2au(0.0);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(0.0);

        // std::size_t n_kx = 5000;
        // std::size_t n_ky = 12;
        // std::size_t n = 61;

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
        // double end = 1.5;
        // Eigen::VectorXd vec = Eigen::VectorXd::LinSpaced(n, start, end);

        // std::ofstream output_file("data/energy.dat");

        // for (auto i = 0; i < n; ++i)
        // {
        //     // sys.mu = meV2au(vec(i));
        //     sys.Bx = T2au(vec(i));
        //     // sys.Bz = T2au(vec(i));

        //     auto evals = calc.eigenvals_sparse_discrete(n_kx, n_ky, 4);
        //     output_file << vec(i) << " " << evals.transpose() / meV2au(1) << std::endl;
        // }
    }

    // discrete system prob den
    {
        // sys.mu = meV2au(4.0);

        // sys.Bx = T2au(1.3);

        // sys.By = T2au(0.0);
        // sys.Bz = T2au(0.0);

        // std::size_t n_kx = 10000;
        // std::size_t n_ky = 12;

        // double energy = meV2au(0.0);

        // printer.printProbDen_sparse_discrete("data/prob_den.dat", n_kx, n_ky, energy);
    }

    // delta vs momentum
    {
        // sys.mu = meV2au(0.0); // doesnt matter for abs delta
        // sys.Bx = T2au(0.0);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(1.0);

        // std::size_t n_k = 201;
        // double k_max = 0.05;

        // Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_k, -k_max, k_max);
        // Eigen::VectorXd ky_vec = Eigen::VectorXd::LinSpaced(n_k, -k_max, k_max);

        // printer.printAbsDelta("data/abs_delta.dat", kx_vec, ky_vec);
        // // sys.printDeltaFromUnitaryTransformation("data/delta.dat", "data/DT.dat", kx_vec, ky_vec);
    }

    // abs delta at given k vs Bz
    {
        // sys.mu = meV2au(0.0); // doesnt matter
        // sys.Bx = T2au(0.0);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(0.0);

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
        //     sys.Bz = T2au(B);
        //     output_file << B << " " << calc.AbsDelta(k.x(), k.y()).transpose() / meV2au(1) << std::endl;
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
        // sys.mu = meV2au(LAOSTO_low_band);

        // sys.Bx = T2au(0.0);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(1.0);

        // printer.printWilsonLoopSpectrum("data/WL.dat", 201, 4001, 0.04);
    }

    // 1D prob_den - two stastes
    {
        // sys.mu = meV2au(LAOSTO_bottom_band);
        // sys.Bx = T2au(0.0);
        // sys.By = T2au(0.0);
        // sys.Bz = T2au(1.0);

        // std::size_t n_ky = 1000;
        // // std::size_t n_ky = 5000;
        // // std::size_t n_ky = 10000;
        // // std::size_t n_ky = 20000;

        // double energy = meV2au(0.0);

        // double kx = 0.005;

        // // std::string filename = "data/" + std::to_string(n_ky) + "C_low_prob_den.dat";
        // std::string filename = "data/" + std::to_string(n_ky) + "C_bot_prob_den.dat";
        // // std::string filename = "data/" + std::to_string(n_ky) + "C_top_prob_den.dat";

        // std::ofstream output_file(filename);

        // // auto [evals, evecs] = calc.eigen_sparse_discrete_ky(kx, n_ky, 2, energy);
        // auto [evals_full, evecs_full] = calc.eigen_discrete_ky(kx, n_ky);
        // auto evals = evals_full.segment(n_ky * 2 * sys.n_bands - 1, 2);
        // auto evecs = evecs_full.block(0, n_ky * 2 * sys.n_bands - 1, 2 * sys.n_bands * n_ky, 2);
        // if (evals(0) > evals(1))
        // {
        //     evals.reverseInPlace();
        //     evecs.col(0).swap(evecs.col(1));
        // }
        // auto prob_den_orbitals = evecs.cwiseAbs2();
        // output_file << "# Eigenvals (meV): " << evals.transpose() / meV2au(1) << " kx = " << kx << " n_ky = " << n_ky << std::endl;
        // for (auto j = 0; j < n_ky; ++j)
        // {
        //     output_file << j << " " << prob_den_orbitals.block(j * 2 * sys.n_bands, 0, 2 * sys.n_bands, 2).colwise().sum() << std::endl;
        // }
    }

    return 0;
}