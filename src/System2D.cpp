#include "System2D.h"
#include "utils.h"

#include <iostream>
#include <fstream>
#include <thread>

System2D::System2D(const H2D &H, const Parameters &p) : _H(H), _p(p)
{
    _n_bands = _H({0, 0}, _p).rows();
}

Eigen::VectorXd System2D::eigenvals(const Point2D &k)
{
    _SAES.compute(_H(k, _p), Eigen::EigenvaluesOnly);
    return _SAES.eigenvalues();
};

Eigen::MatrixXcd System2D::eigenvecs(const Point2D &k)
{
    _SAES.compute(_H(k, _p));
    return _SAES.eigenvectors();
};

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> System2D::eigen(const Point2D &k)
{
    _SAES.compute(_H(k, _p));
    return std::make_tuple(_SAES.eigenvalues(), _SAES.eigenvectors());
};

void System2D::printBandStructure(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
    {
        for (auto ky : ky_vec)
        {
            _SAES.compute(_H({kx, ky}, _p), Eigen::EigenvaluesOnly);
            output_file << kx << " " << ky << " " << _SAES.eigenvalues().transpose() / meV2au(1.0) << "\n";
        }
    }
}

void System2D::printBandStructureSlice(const std::string &output_filename, const Eigen::VectorXd &k_vec, int axis, double k0)
{
    std::ofstream output_file(output_filename);

    if (axis == 0)
    {
        for (auto kx : k_vec)
        {
            _SAES.compute(_H({kx, k0}, _p), Eigen::EigenvaluesOnly);
            output_file << kx << " " << _SAES.eigenvalues().transpose() / meV2au(1.0) << "\n";
        }
    }
    else if (axis == 1)
    {
        for (auto ky : k_vec)
        {
            _SAES.compute(_H({k0, ky}, _p), Eigen::EigenvaluesOnly);
            output_file << ky << " " << _SAES.eigenvalues().transpose() / meV2au(1.0) << "\n";
        }
    }
    else
    {
        std::cerr << "Invalid axis. Choose 0 or 1." << std::endl;
    }
}

void System2D::printGap(const std::string &output_filename, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec)
{
    std::ofstream output_file(output_filename);

    for (auto kx : kx_vec)
    {
        for (auto ky : ky_vec)
        {
            _SAES.compute(_H({kx, ky}, _p), Eigen::EigenvaluesOnly);
            output_file << kx << " " << ky << " " << (_SAES.eigenvalues().tail(_n_bands / 2) - _SAES.eigenvalues().head(_n_bands / 2).reverse()).transpose() / meV2au(1.0) << "\n";
        }
    }
}

void System2D::printGapAlongContour(const std::string &output_filename, const std::vector<Point2D> &contour)
{
    std::ofstream output_file(output_filename);

    for (auto k : contour)
    {
        _SAES.compute(_H(k, _p), Eigen::EigenvaluesOnly);
        output_file << k.x() << " " << k.y() << " " << (_SAES.eigenvalues().tail(_n_bands / 2) - _SAES.eigenvalues().head(_n_bands / 2).reverse()).transpose() / meV2au(1.0) << "\n";
    }
}

std::vector<std::vector<Point2D>> System2D::findFSContours(double E, double dk, double eps, const Point2D kx_range, std::size_t n_kx)
{
    std::vector<std::vector<Point2D>> FS_contours;

    // search for kx such that E is an eigenvalue of H(kx, 0)
    Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_kx, kx_range(0), kx_range(1));

    _SAES.compute(_H({kx_vec(0), 0.0}, _p), Eigen::EigenvaluesOnly);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> SAES_next;

    std::vector<std::size_t> kx_idx; // indices of kx_vec where the contour is found
    std::vector<std::size_t> E_idx;  // indices of bands creating the contour

    for (auto ie = 0; ie < _n_bands; ++ie) // contour index
    {
        for (auto i = (kx_idx.empty() ? 0 : kx_idx.back()); i < kx_vec.size() - 1; ++i)
        {
            switch (i % 2)
            {
            case 0:
                SAES_next.compute(_H({kx_vec(i + 1), 0.0}, _p), Eigen::EigenvaluesOnly);
                break;
            case 1:
                _SAES.compute(_H({kx_vec(i + 1), 0.0}, _p), Eigen::EigenvaluesOnly);
                break;
            }

            if ((_SAES.eigenvalues()(ie) - E) * (SAES_next.eigenvalues()(ie) - E) < 0)
            {
                // std::cout << "contour " << ic << " found at kx = " << kx_vec(i) << std::endl;
                kx_idx.push_back(i);
                E_idx.push_back(ie);
                break;
            }
        }
    }

    auto R_mat = [](double phi) -> Eigen::Matrix2d
    {
        Eigen::Matrix2d R;
        R << std::cos(phi), -std::sin(phi), std::sin(phi), std::cos(phi);
        return R;
    };

    // begin looking for contours in -ky direction
    for (auto ic = 0; ic < kx_idx.size(); ++ic)
    {
        // find exact start kx with regula falsi
        double kxl = kx_vec[kx_idx[ic]];
        double kxr = kx_vec[kx_idx[ic] + 1];

        double El = _SAES.compute(_H({kxl, 0.0}, _p), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);
        double Er = _SAES.compute(_H({kxr, 0.0}, _p), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);

        double kx_exact = 0.0;
        double E_exact = 0.0;

        do
        {
            kx_exact = (El * kxr - Er * kxl) / (El - Er);
            E_exact = _SAES.compute(_H({kx_exact, 0.0}, _p), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);

            // if (std::abs(E_exact - E) < eps)
            //     break;

            if (El * E_exact < 0)
            {
                kxr = kx_exact;
                Er = E_exact;
            }
            else
            {
                kxl = kx_exact;
                El = E_exact;
            }

        } while (std::abs(E_exact - E) > eps);

        // std::cout << "contour " << ic << " found at kx_exact = " << kx_exact << std::endl;

        // first point of the contour
        std::vector<Point2D> contour;
        contour.push_back(Point2D{kx_exact, 0.0});
        contour.reserve(1e5); // try to avoid reallocations

        // find next points of the contour
        double phil = 0.0;
        double phir = 0.0;

        Eigen::Matrix2d Rl;
        Eigen::Matrix2d Rr;

        Point2D k_diff;
        Point2D kl;
        Point2D kr;

        double phi_next = 0.0;
        Eigen::Matrix2d R_next;
        Point2D k_next = {0.0, 0.0};
        double E_next = 0.0;

        do
        {
            phil = M_PI_4;
            phir = -M_PI_4;

            Rl = R_mat(phil);
            Rr = Rl;
            Rr(0, 1) *= -1;
            Rr(1, 0) *= -1;

            k_diff = contour.size() > 1 ? contour.end()[-1] - contour.end()[-2] : Point2D{0.0, -dk};

            kl = contour.end()[-1] + Rl * k_diff;
            kr = contour.end()[-1] + Rr * k_diff;

            El = _SAES.compute(_H(kl, _p), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);
            Er = _SAES.compute(_H(kr, _p), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);

            // find next k - regula falsi with phi
            do
            {
                phi_next = (El * phir - Er * phil) / (El - Er);
                R_next = R_mat(phi_next);
                k_next = contour.end()[-1] + R_next * k_diff;
                E_next = _SAES.compute(_H(k_next, _p), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);

                // if (std::abs(E_next - E) < eps)
                //     break;

                if (El * E_exact < 0)
                {
                    phir = phi_next;
                    Er = E_next;
                }
                else
                {
                    phil = phi_next;
                    El = E_next;
                }

            } while (std::abs(E_next - E) > eps);

            contour.push_back(k_next);
        } while (contour.size() < 3 || (contour.end()[-1] - contour[0]).norm() > dk);

        // close the contour
        if (contour.back().y() < 0)
            contour.back() = contour[0];
        else
            contour.push_back(contour[0]);

        FS_contours.push_back(contour);
    }

    return FS_contours;
}

void System2D::printBerryCurvature(const std::string &output_filename, const H2D &dvxH, const H2D &dvyH, const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec)
{
    std::ofstream output_file(output_filename);

    Eigen::VectorXd BC(_n_bands); // Berry curvature

    Eigen::MatrixXcd dvxH_mat(_n_bands, _n_bands);
    Eigen::MatrixXcd dvyH_mat(_n_bands, _n_bands);

    Eigen::dcomplex v1x, v1y, v2x, v2y;

    for (auto kx : kx_vec)
    {
        for (auto ky : ky_vec)
        {
            _SAES.compute(_H({kx, ky}, _p));

            dvxH_mat = dvxH({kx, ky}, _p);
            dvyH_mat = dvyH({kx, ky}, _p);

            BC.setZero();

            for (auto bi = 0; bi < _n_bands; ++bi) // bi - band index
            {

                for (auto bj = 0; bj < _n_bands; ++bj)
                {
                    if (bj == bi)
                        continue;

                    v1x = _SAES.eigenvectors().col(bi).dot(dvxH_mat * _SAES.eigenvectors().col(bj));
                    v1y = _SAES.eigenvectors().col(bi).dot(dvyH_mat * _SAES.eigenvectors().col(bj));

                    v2x = _SAES.eigenvectors().col(bj).dot(dvxH_mat * _SAES.eigenvectors().col(bi));
                    v2y = _SAES.eigenvectors().col(bj).dot(dvyH_mat * _SAES.eigenvectors().col(bi));

                    BC(bi) -= std::imag(v1x * v2y - v2x * v1y) / (_SAES.eigenvalues()(bi) - _SAES.eigenvalues()(bj)) / (_SAES.eigenvalues()(bi) - _SAES.eigenvalues()(bj));
                }
            }

            output_file << kx << " " << ky << " " << BC.transpose() * 1e6 << "\n";
        }
    }
}

// Multithreaded version - second try
Eigen::VectorXd System2D::calcChernNumbers(const Eigen::Vector2<std::size_t> &n, double kmax)
{
    if (n.x() % 2 != 0 || n.y() % 2 != 0)
    {
        std::cerr << "Number of k-points must be even in each direction. Returning empty vector." << std::endl;
        return Eigen::VectorXd();
    }

    // BZ mesh
    Eigen::VectorXd kx;
    Eigen::VectorXd ky;

    if (kmax >= M_PI)
    {
        kx = Eigen::VectorXd::LinSpaced(n.x(), -M_PI, M_PI);
        ky = Eigen::VectorXd::LinSpaced(n.y(), -M_PI, M_PI);

        // account for periodic boundary conditions
        kx(n.x() - 1) = kx(0);
        ky(n.y() - 1) = ky(0);
    }
    else
    {
        kx = Eigen::VectorXd::LinSpaced(n.x(), -kmax, kmax);
        ky = Eigen::VectorXd::LinSpaced(n.y(), -kmax, kmax);
    }

    return calcChernNumbersWithCustomGrid(kx, ky);
}

Eigen::VectorXd System2D::calcChernNumbersDenserCenter(std::size_t n_dense, std::size_t n_sparse, double k_val)
{
    // dense region
    Eigen::VectorXd kc = Eigen::VectorXd::LinSpaced(n_dense, -k_val, k_val);

    // without the dense region
    Eigen::VectorXd ko = Eigen::VectorXd::LinSpaced(n_sparse, -M_PI, M_PI);
    ko(ko.size() - 1) = -M_PI;

    // before and after the dense region
    std::size_t n = ((1.0 - k_val / M_PI) * n_sparse) / 2.0;
    if (n % 2 != 0)
        n++;

    Eigen::VectorXd k_before = Eigen::VectorXd::LinSpaced(n, -M_PI, -k_val);
    Eigen::VectorXd k_after = Eigen::VectorXd::LinSpaced(n, k_val, M_PI);
    k_after(k_after.size() - 1) = -M_PI;

    // kx next to the dense region
    std::size_t m = k_val / M_PI * n_sparse;
    if (m % 2 != 0)
        m++;
    Eigen::VectorXd kxn = Eigen::VectorXd::LinSpaced(m, -k_val, k_val);

    // std::cout << "n = " << n << ", m = " << m << std::endl;
    // std::cout << "dense region: " << kc(0) << "," << kc(0) << " - " << kc(kc.size() - 1) << "," << kc(kc.size() - 1) << std::endl;
    // std::cout << "top bar: " << k_before(0) << "," << ko(0) << " - " << k_before(k_before.size() - 1) << "," << ko(ko.size() - 1) << std::endl;
    // std::cout << "bottom bar: " << k_after(0) << "," << ko(0) << " - " << k_after(k_after.size() - 1) << "," << ko(ko.size() - 1) << std::endl;
    // std::cout << "left smaller bar: " << kxn(0) << "," << k_before(0) << " - " << kxn(kxn.size() - 1) << "," << k_before(k_before.size() - 1) << std::endl;
    // std::cout << "right smaller bar: " << kxn(0) << "," << k_after(0) << " - " << kxn(kxn.size() - 1) << "," << k_after(k_after.size() - 1) << std::endl;

    return calcChernNumbersWithCustomGrid(kc, kc)          // dense region
           + calcChernNumbersWithCustomGrid(k_before, ko)  // top bar
           + calcChernNumbersWithCustomGrid(k_after, ko)   // botom bar
           + calcChernNumbersWithCustomGrid(kxn, k_before) // left smaller bar
           + calcChernNumbersWithCustomGrid(kxn, k_after); // right smaller bar
}

Eigen::VectorXd System2D::calcChernNumbersWithCustomGrid(const Eigen::VectorXd &kx, const Eigen::VectorXd &ky)
{
    std::size_t nx = kx.size();
    std::size_t ny = ky.size();

    if (nx % 2 != 0 || ny % 2 != 0)
    {
        std::cerr << "Number of k-points must be even in each direction. Returning empty vector." << std::endl;
        return Eigen::VectorXd();
    }

    // BZ mesh
    Eigen::VectorXd kx_t = kx.head(nx / 2 + 1);
    Eigen::VectorXd kx_b = kx.tail(nx / 2);

    Eigen::VectorXd ky_l = ky.head(ny / 2 + 1);
    Eigen::VectorXd ky_r = ky.tail(ny / 2);

    // with rows storing - original
    // auto processGridSection = [this](const Eigen::VectorXd &kx, const Eigen::VectorXd &ky, Eigen::VectorXd &CN) -> void
    // {
    //     CN.setZero();

    //     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> SAES;

    //     std::size_t rows = kx.size();
    //     std::size_t cols = ky.size();

    //     Eigen::VectorX<Eigen::MatrixXcd> cr_evecs(ky.size()); // Current row eigenvectors

    //     // Compute eigenvectors for the first row
    //     for (auto j = 0; j < cols; ++j)
    //     {
    //         SAES.compute(_H({kx(0), ky(j)}, _p));
    //         cr_evecs(j) = SAES.eigenvectors();
    //     }

    //     Eigen::MatrixXcd ni_j_evecs(_n_bands, _n_bands); // Next i, j eigenvectors

    //     Eigen::VectorXd BCdkxdky(_n_bands); // Berry curvature * dkx * dky
    //     BCdkxdky.setZero();

    //     Eigen::dcomplex yp, xp, ym, xm; // eigenvectors dot products, yp in y+ direction, xm in x- direction etc.

    //     std::size_t ni = 1; // Next i
    //     std::size_t nj = 1; // Next j

    //     for (auto i = 0; i < rows - 1; ++i)
    //     {
    //         ni = i + 1;

    //         SAES.compute(_H({kx(ni), ky(0)}, _p));
    //         ni_j_evecs = SAES.eigenvectors();

    //         for (auto j = 0; j < cols - 1; ++j)
    //         {
    //             nj = j + 1;

    //             SAES.compute(_H({kx(ni), ky(nj)}, _p)); // contains eigenvectors and eigenvalues for next i, next j

    //             for (auto bi = 0; bi < _n_bands; ++bi) // bi - band index
    //             {
    //                 yp = cr_evecs(j).col(bi).dot(cr_evecs(nj).col(bi));         // psi_{i,j} dot psi_{i,j+1}
    //                 xp = cr_evecs(nj).col(bi).dot(SAES.eigenvectors().col(bi)); // psi_{i,j+1} dot psi_{i+1,j+1}
    //                 ym = SAES.eigenvectors().col(bi).dot(ni_j_evecs.col(bi));   // psi_{i+1,j+1} dot psi_{i+1,j}
    //                 xm = ni_j_evecs.col(bi).dot(cr_evecs(j).col(bi));           // psi_{i+1,j} dot psi_{i,j}

    //                 BCdkxdky(bi) = -std::arg(yp * xp * ym * xm); // std::arg uses std::atan2, maybe negative?
    //             }

    //             CN += BCdkxdky;

    //             cr_evecs(j) = std::move(ni_j_evecs); // set used current row i,j eigenvectors to ni eigenvectors for next i
    //             ni_j_evecs = SAES.eigenvectors();    // set ni,j eigenvectors to ni,nj eigenvectors for next j
    //         }

    //         cr_evecs(cols - 1) = std::move(ni_j_evecs); // account for last in row
    //     }

    //     CN /= (2.0 * M_PI);
    // };

    // with rows storing and jumping in rows
    auto processGridSection = [this](const Eigen::VectorXd &kx, const Eigen::VectorXd &ky, Eigen::VectorXd &CN) -> void
    {
        CN.setZero();

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> SAES;

        std::size_t rows = kx.size();
        std::size_t cols = ky.size();

        Eigen::VectorX<Eigen::MatrixXcd> i_evecs(cols);  // Current row eigenvectors
        Eigen::VectorX<Eigen::MatrixXcd> ni_evecs(cols); // Next row eigenvectors

        Eigen::VectorX<Eigen::MatrixXcd> *i_evecs_ptr = &i_evecs;
        Eigen::VectorX<Eigen::MatrixXcd> *ni_evecs_ptr = &ni_evecs;

        // Compute eigenvectors for the first row
        for (auto j = 0; j < cols; ++j)
        {
            i_evecs(j) = SAES.compute(_H({kx(0), ky(j)}, _p)).eigenvectors();
        }

        Eigen::VectorXd BCdkxdky(_n_bands); // Berry curvature * dkx * dky
        BCdkxdky.setZero();

        Eigen::dcomplex yp, xp, ym, xm; // eigenvectors dot products, yp in y+ direction, xm in x- direction etc.

        for (auto i = 0; i < rows - 1; ++i)
        {
            (*ni_evecs_ptr)(0) = SAES.compute(_H({kx(i + 1), ky(0)}, _p)).eigenvectors();

            for (auto j = 0; j < cols - 1; ++j)
            {
                (*ni_evecs_ptr)(j + 1) = SAES.compute(_H({kx(i + 1), ky(j + 1)}, _p)).eigenvectors(); // contains eigenvectors and eigenvalues for next i, next j

                for (auto bi = 0; bi < _n_bands; ++bi) // bi - band index
                {
                    yp = (*i_evecs_ptr)(j).col(bi).dot((*i_evecs_ptr)(j + 1).col(bi));      // psi_{i,j} dot psi_{i,j+1}
                    xp = (*i_evecs_ptr)(j + 1).col(bi).dot((*ni_evecs_ptr)(j + 1).col(bi)); // psi_{i,j+1} dot psi_{i+1,j+1}
                    ym = (*ni_evecs_ptr)(j + 1).col(bi).dot((*ni_evecs_ptr)(j).col(bi));    // psi_{i+1,j+1} dot psi_{i+1,j}
                    xm = (*ni_evecs_ptr)(j).col(bi).dot((*i_evecs_ptr)(j).col(bi));         // psi_{i+1,j} dot psi_{i,j}

                    BCdkxdky(bi) = -std::arg(yp * xp * ym * xm); // std::arg uses std::atan2, maybe negative?
                    // BCdkxdky(bi) = -std::log(yp * xp * ym * xm).imag(); // std::arg uses std::atan2, maybe negative?
                }

                CN += BCdkxdky;
            }

            std::swap(i_evecs_ptr, ni_evecs_ptr);
        }

        CN /= (2.0 * M_PI);
    };

    // // simpler version without rows storing
    // auto processGridSection = [this](const Eigen::VectorXd &kx, const Eigen::VectorXd &ky, Eigen::VectorXd &CN) -> void
    // {
    //     CN.setZero();

    //     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> SAES_i_j;
    //     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> SAES_ni_j;
    //     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> SAES_i_nj;
    //     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> SAES_ni_nj;

    //     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> *SAES_i_j_ptr = &SAES_i_j;
    //     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> *SAES_ni_j_ptr = &SAES_ni_j;
    //     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> *SAES_i_nj_ptr = &SAES_i_nj;
    //     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> *SAES_ni_nj_ptr = &SAES_ni_nj;

    //     std::size_t rows = kx.size();
    //     std::size_t cols = ky.size();

    //     Eigen::VectorXd BCdkxdky(_n_bands); // Berry curvature * dkx * dky
    //     BCdkxdky.setZero();

    //     Eigen::dcomplex yp, xp, ym, xm; // eigenvectors dot products, yp in y+ direction, xm in x- direction etc.

    //     for (auto i = 0; i < rows - 1; ++i)
    //     {
    //         SAES_i_j.compute(_H({kx(i), ky(0)}, _p));
    //         SAES_ni_j.compute(_H({kx(i + 1), ky(0)}, _p));

    //         SAES_i_j_ptr = &SAES_i_j;
    //         SAES_ni_j_ptr = &SAES_ni_j;
    //         SAES_i_nj_ptr = &SAES_i_nj;
    //         SAES_ni_nj_ptr = &SAES_ni_nj;

    //         for (auto j = 0; j < cols - 1; ++j)
    //         {
    //             SAES_i_nj_ptr->compute(_H({kx(i), ky(j + 1)}, _p));
    //             SAES_ni_nj_ptr->compute(_H({kx(i + 1), ky(j + 1)}, _p));

    //             // #pragma omp parallel for reduction(+ : BCdkxdky)
    //             for (auto bi = 0; bi < _n_bands; ++bi) // bi - band index
    //             {
    //                 yp = SAES_i_j_ptr->eigenvectors().col(bi).dot(SAES_i_nj_ptr->eigenvectors().col(bi));   // psi_{i,j} dot psi_{i,j+1}
    //                 xp = SAES_i_nj_ptr->eigenvectors().col(bi).dot(SAES_ni_nj_ptr->eigenvectors().col(bi)); // psi_{i,j+1} dot psi_{i+1,j+1}
    //                 ym = SAES_ni_nj_ptr->eigenvectors().col(bi).dot(SAES_ni_j_ptr->eigenvectors().col(bi)); // psi_{i+1,j+1} dot psi_{i+1,j}
    //                 xm = SAES_ni_j_ptr->eigenvectors().col(bi).dot(SAES_i_j_ptr->eigenvectors().col(bi));   // psi_{i+1,j} dot psi_{i,j}

    //                 BCdkxdky(bi) = -std::arg(yp * xp * ym * xm); // std::arg uses std::atan2, maybe negative?
    //             }

    //             CN += BCdkxdky;

    //             std::swap(SAES_i_j_ptr, SAES_i_nj_ptr);
    //             std::swap(SAES_ni_j_ptr, SAES_ni_nj_ptr);
    //         }
    //     }

    //     CN /= (2.0 * M_PI);
    // };

    Eigen::VectorXd CN_tl(_n_bands);
    Eigen::VectorXd CN_tr(_n_bands);
    Eigen::VectorXd CN_bl(_n_bands);
    Eigen::VectorXd CN_br(_n_bands);

    std::thread tl(processGridSection, kx_t, ky_l, std::ref(CN_tl));
    std::thread tr(processGridSection, kx_t, ky_r, std::ref(CN_tr));
    std::thread bl(processGridSection, kx_b, ky_l, std::ref(CN_bl));
    // std::thread br(processGridSection, kx_b, ky_r, std::ref(CN_br));

    // bottom right using main thread
    processGridSection(kx_b, ky_r, CN_br);

    tl.join();
    tr.join();
    bl.join();
    // br.join();

    return CN_tl + CN_tr + CN_bl + CN_br;
}

// w zasadzie dziala, ale jest wolne i dzwoni
Eigen::VectorXd System2D::calcChernNumbersFromBC(const Eigen::Vector2<std::size_t> &n, const H2D &dvxH, const H2D &dvyH)
{

    if (n.x() % 2 != 0 || n.y() % 2 != 0)
    {
        std::cerr << "Number of k-points must be even in each direction. Returning empty vector." << std::endl;
        return Eigen::VectorXd();
    }

    Eigen::VectorXd kx = Eigen::VectorXd::LinSpaced(n.x(), -M_PI, M_PI);
    Eigen::VectorXd ky = Eigen::VectorXd::LinSpaced(n.y(), -M_PI, M_PI);

    kx = kx.head(n.x() - 1).eval();
    ky = ky.head(n.y() - 1).eval();

    Eigen::VectorXd kx_t = kx.head(n.x() / 2 + 1);
    Eigen::VectorXd kx_b = kx.tail(n.x() / 2);

    Eigen::VectorXd ky_l = ky.head(n.y() / 2 + 1);
    Eigen::VectorXd ky_r = ky.tail(n.y() / 2);

    auto processGridSection = [this, dvxH, dvyH](const Eigen::VectorXd &kx_vec, const Eigen::VectorXd &ky_vec, Eigen::VectorXd &CN) -> void
    {
        CN.setZero();

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> SAES;

        Eigen::VectorXd BC(_n_bands); // Berry curvature

        Eigen::MatrixXcd dvxH_mat(_n_bands, _n_bands);
        Eigen::MatrixXcd dvyH_mat(_n_bands, _n_bands);

        Eigen::dcomplex v1x, v1y, v2x, v2y;

        for (auto kx : kx_vec)
        {
            for (auto ky : ky_vec)
            {
                SAES.compute(_H({kx, ky}, _p));

                dvxH_mat = dvxH({kx, ky}, _p);
                dvyH_mat = dvyH({kx, ky}, _p);

                BC.setZero();

                for (auto bi = 0; bi < _n_bands; ++bi) // bi - band index
                {

                    for (auto bj = 0; bj < _n_bands; ++bj)
                    {
                        if (bj == bi)
                            continue;

                        v1x = SAES.eigenvectors().col(bi).dot(dvxH_mat * SAES.eigenvectors().col(bj));
                        v1y = SAES.eigenvectors().col(bi).dot(dvyH_mat * SAES.eigenvectors().col(bj));

                        v2x = SAES.eigenvectors().col(bj).dot(dvxH_mat * SAES.eigenvectors().col(bi));
                        v2y = SAES.eigenvectors().col(bj).dot(dvyH_mat * SAES.eigenvectors().col(bi));

                        BC(bi) -= std::imag(v1x * v2y - v2x * v1y) / (SAES.eigenvalues()(bi) - SAES.eigenvalues()(bj)) / (SAES.eigenvalues()(bi) - SAES.eigenvalues()(bj));
                    }
                }

                CN += BC;
            }
        }

        CN *= (kx_vec(1) - kx_vec(0)) * (ky_vec(1) - ky_vec(0)) / (2.0 * M_PI);
    };

    Eigen::VectorXd CN_tl(_n_bands);
    Eigen::VectorXd CN_tr(_n_bands);
    Eigen::VectorXd CN_bl(_n_bands);
    Eigen::VectorXd CN_br(_n_bands);

    std::thread tl(processGridSection, kx_t, ky_l, std::ref(CN_tl));
    std::thread tr(processGridSection, kx_t, ky_r, std::ref(CN_tr));
    std::thread bl(processGridSection, kx_b, ky_l, std::ref(CN_bl));
    // std::thread br(processGridSection, kx_b, ky_r, std::ref(CN_br));

    // bottom right using main thread
    processGridSection(kx_b, ky_r, CN_br);

    tl.join();
    tr.join();
    bl.join();
    // br.join();

    return CN_tl + CN_tr + CN_bl + CN_br;
}