#include "System2DCalculations.h"
#include "utils.h"

#include <iostream>
#include <fstream>
#include <thread>

Hamiltonian System2DCalculations::dHdkx(double kx, double ky, double dk) const
{
    return (_sys.HBdG(kx + dk, ky) - _sys.HBdG(kx - dk, ky)) / (2.0 * dk);
}

Hamiltonian System2DCalculations::dHdky(double kx, double ky, double dk) const
{
    return (_sys.HBdG(kx, ky + dk) - _sys.HBdG(kx, ky - dk)) / (2.0 * dk);
}

Hamiltonian System2DCalculations::dHdkx_normal(double kx, double ky, double dk) const
{
    return (_sys.Hk(kx + dk, ky) - _sys.Hk(kx - dk, ky)) / (2.0 * dk);
}

Hamiltonian System2DCalculations::dHdky_normal(double kx, double ky, double dk) const
{
    return (_sys.Hk(kx, ky + dk) - _sys.Hk(kx, ky - dk)) / (2.0 * dk);
}

Eigen::VectorXd System2DCalculations::generate_k_vec(std::size_t n_dense, std::size_t n_sparse, double k_val) const
{

    Eigen::VectorXd k_vec(n_dense + 2 * n_sparse + 1);
    k_vec.segment(0, n_sparse + 1) = Eigen::VectorXd::LinSpaced(n_sparse + 1, -M_PI, -k_val);
    k_vec.segment(n_sparse, n_dense + 1) = Eigen::VectorXd::LinSpaced(n_dense + 1, -k_val, k_val);
    k_vec.segment(n_sparse + n_dense, n_sparse + 1) = Eigen::VectorXd::LinSpaced(n_sparse + 1, k_val, M_PI);
    k_vec(k_vec.size() - 1) = -M_PI;

    return k_vec;
}

Eigen::VectorXd System2DCalculations::eigenvals(double kx, double ky)
{
    return _SAES.compute(_sys.HBdG(kx, ky), Eigen::EigenvaluesOnly).eigenvalues();
};

Eigen::VectorXd System2DCalculations::eigenvals_normal(double kx, double ky)
{
    return _SAES.compute(_sys.Hk(kx, ky), Eigen::EigenvaluesOnly).eigenvalues();
};

Eigen::VectorXd System2DCalculations::eigenvals_discrete_kx(std::size_t n_kx, double ky)
{
    return _SAES.compute(_sys.HBdG_discrete_kx(n_kx, ky), Eigen::EigenvaluesOnly).eigenvalues();
};

Eigen::VectorXd System2DCalculations::eigenvals_discrete_kx_normal(std::size_t n_kx, double ky)
{
    return _SAES.compute(_sys.Hk_discrete_kx(n_kx, ky), Eigen::EigenvaluesOnly).eigenvalues();
};

Eigen::VectorXd System2DCalculations::eigenvals_discrete_ky(double kx, std::size_t n_ky)
{
    return _SAES.compute(_sys.HBdG_discrete_ky(kx, n_ky), Eigen::EigenvaluesOnly).eigenvalues();
};

Eigen::VectorXd System2DCalculations::eigenvals_discrete_ky_normal(double kx, std::size_t n_ky)
{
    return _SAES.compute(_sys.Hk_discrete_ky(kx, n_ky), Eigen::EigenvaluesOnly).eigenvalues();
};

Eigen::VectorXd System2DCalculations::eigenvals_discrete(std::size_t n_kx, std::size_t n_ky)
{
    return _SAES.compute(_sys.HBdG_discrete(n_kx, n_ky), Eigen::EigenvaluesOnly).eigenvalues();
};

Eigen::VectorXd System2DCalculations::eigenvals_discrete_normal(std::size_t n_kx, std::size_t n_ky)
{
    return _SAES.compute(_sys.Hk_discrete(n_kx, n_ky), Eigen::EigenvaluesOnly).eigenvalues();
};

Eigen::VectorXd System2DCalculations::eigenvals_sparse_discrete_kx(std::size_t n_kx, double ky, std::size_t n_eigs, double sigma)
{
    return eigenvals_sparse(_sys.HBdG_discrete_kx_sparse(n_kx, ky), n_eigs, sigma);
};

Eigen::VectorXd System2DCalculations::eigenvals_sparse_discrete_kx_normal(std::size_t n_kx, double ky, std::size_t n_eigs, double sigma)
{
    return eigenvals_sparse(_sys.Hk_discrete_kx_sparse(n_kx, ky), n_eigs, sigma);
};

Eigen::VectorXd System2DCalculations::eigenvals_sparse_discrete_ky(double kx, std::size_t n_ky, std::size_t n_eigs, double sigma)
{
    return eigenvals_sparse(_sys.HBdG_discrete_ky_sparse(kx, n_ky), n_eigs, sigma);
};

Eigen::VectorXd System2DCalculations::eigenvals_sparse_discrete_ky_normal(double kx, std::size_t n_ky, std::size_t n_eigs, double sigma)
{
    return eigenvals_sparse(_sys.Hk_discrete_ky_sparse(kx, n_ky), n_eigs, sigma);
};

Eigen::VectorXd System2DCalculations::eigenvals_sparse_discrete(std::size_t n_kx, std::size_t n_ky, std::size_t n_eigs, double sigma)
{
    return eigenvals_sparse(_sys.HBdG_discrete_sparse(n_kx, n_ky), n_eigs, sigma);
};

Eigen::VectorXd System2DCalculations::eigenvals_sparse_discrete_normal(std::size_t n_kx, std::size_t n_ky, std::size_t n_eigs, double sigma)
{
    return eigenvals_sparse(_sys.Hk_discrete_sparse(n_kx, n_ky), n_eigs, sigma);
};

Eigen::MatrixXcd System2DCalculations::eigenvecs(double kx, double ky)
{
    return _SAES.compute(_sys.HBdG(kx, ky)).eigenvectors();
};

Eigen::MatrixXcd System2DCalculations::eigenvecs_normal(double kx, double ky)
{
    return _SAES.compute(_sys.Hk(kx, ky)).eigenvectors();
};

Eigen::MatrixXcd System2DCalculations::eigenvecs_discrete_kx(std::size_t n_kx, double ky)
{
    return _SAES.compute(_sys.HBdG_discrete_kx(n_kx, ky)).eigenvectors();
};

Eigen::MatrixXcd System2DCalculations::eigenvecs_discrete_kx_normal(std::size_t n_kx, double ky)
{
    return _SAES.compute(_sys.Hk_discrete_kx(n_kx, ky)).eigenvectors();
};

Eigen::MatrixXcd System2DCalculations::eigenvecs_discrete_ky(double kx, std::size_t n_ky)
{
    return _SAES.compute(_sys.HBdG_discrete_ky(kx, n_ky)).eigenvectors();
};

Eigen::MatrixXcd System2DCalculations::eigenvecs_discrete_ky_normal(double kx, std::size_t n_ky)
{
    return _SAES.compute(_sys.Hk_discrete_ky(kx, n_ky)).eigenvectors();
};

Eigen::MatrixXcd System2DCalculations::eigenvecs_discrete(std::size_t n_kx, std::size_t n_ky)
{
    return _SAES.compute(_sys.HBdG_discrete(n_kx, n_ky)).eigenvectors();
};

Eigen::MatrixXcd System2DCalculations::eigenvecs_discrete_normal(std::size_t n_kx, std::size_t n_ky)
{
    return _SAES.compute(_sys.Hk_discrete(n_kx, n_ky)).eigenvectors();
};

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> System2DCalculations::eigen(double kx, double ky)
{
    _SAES.compute(_sys.HBdG(kx, ky));
    return std::make_pair(_SAES.eigenvalues(), _SAES.eigenvectors());
};

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> System2DCalculations::eigen_normal(double kx, double ky)
{
    _SAES.compute(_sys.Hk(kx, ky));
    return std::make_pair(_SAES.eigenvalues(), _SAES.eigenvectors());
};

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> System2DCalculations::eigen_discrete_kx(std::size_t n_kx, double ky)
{
    _SAES.compute(_sys.HBdG_discrete_kx(n_kx, ky));
    return std::make_pair(_SAES.eigenvalues(), _SAES.eigenvectors());
};

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> System2DCalculations::eigen_discrete_kx_normal(std::size_t n_kx, double ky)
{
    _SAES.compute(_sys.Hk_discrete_kx(n_kx, ky));
    return std::make_pair(_SAES.eigenvalues(), _SAES.eigenvectors());
};

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> System2DCalculations::eigen_discrete_ky(double kx, std::size_t n_ky)
{
    _SAES.compute(_sys.HBdG_discrete_ky(kx, n_ky));
    return std::make_pair(_SAES.eigenvalues(), _SAES.eigenvectors());
};

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> System2DCalculations::eigen_discrete_ky_normal(double kx, std::size_t n_ky)
{
    _SAES.compute(_sys.Hk_discrete_ky(kx, n_ky));
    return std::make_pair(_SAES.eigenvalues(), _SAES.eigenvectors());
};

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> System2DCalculations::eigen_discrete(std::size_t n_kx, std::size_t n_ky)
{
    _SAES.compute(_sys.HBdG_discrete(n_kx, n_ky));
    return std::make_pair(_SAES.eigenvalues(), _SAES.eigenvectors());
};

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> System2DCalculations::eigen_discrete_normal(std::size_t n_kx, std::size_t n_ky)
{
    _SAES.compute(_sys.Hk_discrete(n_kx, n_ky));
    return std::make_pair(_SAES.eigenvalues(), _SAES.eigenvectors());
};

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> System2DCalculations::eigen_sparse_discrete_kx(std::size_t n_kx, double ky, std::size_t n_eigs, double sigma)
{
    return eigen_sparse(_sys.HBdG_discrete_kx_sparse(n_kx, ky), n_eigs, sigma);
};

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> System2DCalculations::eigen_sparse_discrete_kx_normal(std::size_t n_kx, double ky, std::size_t n_eigs, double sigma)
{
    return eigen_sparse(_sys.Hk_discrete_kx_sparse(n_kx, ky), n_eigs, sigma);
};

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> System2DCalculations::eigen_sparse_discrete_ky(double kx, std::size_t n_ky, std::size_t n_eigs, double sigma)
{
    return eigen_sparse(_sys.HBdG_discrete_ky_sparse(kx, n_ky), n_eigs, sigma);
}

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> System2DCalculations::eigen_sparse_discrete_ky_normal(double kx, std::size_t n_ky, std::size_t n_eigs, double sigma)
{
    return eigen_sparse(_sys.Hk_discrete_ky_sparse(kx, n_ky), n_eigs, sigma);
}

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> System2DCalculations::eigen_sparse_discrete(std::size_t n_kx, std::size_t n_ky, std::size_t n_eigs, double sigma)
{
    return eigen_sparse(_sys.HBdG_discrete_sparse(n_kx, n_ky), n_eigs, sigma);
}

std::pair<Eigen::VectorXd, Eigen::MatrixXcd> System2DCalculations::eigen_sparse_discrete_normal(std::size_t n_kx, std::size_t n_ky, std::size_t n_eigs, double sigma)
{
    return eigen_sparse(_sys.Hk_discrete_sparse(n_kx, n_ky), n_eigs, sigma);
}

Eigen::VectorXd System2DCalculations::AbsDelta(double kx, double ky)
{
    Eigen::VectorXd abs_deltas(_sys.n_bands);

    double mu_original = _sys.mu;
    _sys.mu = 0.0;

    Eigen::VectorXd evals_mu_zero = eigenvals_normal(kx, ky);

    for (auto i = 0; i < _sys.n_bands; ++i)
    {
        _sys.mu = evals_mu_zero(i);
        Eigen::VectorXd evals = eigenvals(kx, ky);
        abs_deltas(i) = (evals(_sys.n_bands) - evals(_sys.n_bands - 1)) / 2.0;
    }

    _sys.mu = mu_original;

    return abs_deltas;
}

// abelowe-dla oddzielonych od siebie pasm
Eigen::VectorXd System2DCalculations::AbelianBerryCurvature(double kx, double ky)
{
    Eigen::VectorXd BC = Eigen::VectorXd::Zero(2 * _sys.n_bands);

    Eigen::MatrixXcd dvxH_mat = dHdkx(kx, ky);
    Eigen::MatrixXcd dvyH_mat = dHdky(kx, ky);

    Eigen::dcomplex v1x, v1y, v2x, v2y;

    // local solver to allow multi-threading
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> SAES(_sys.HBdG(kx, ky));

    for (auto bi = 0; bi < 2 * _sys.n_bands; ++bi) // bi - band index
    {
        for (auto bj = 0; bj < 2 * _sys.n_bands; ++bj)
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

    return BC;
}

double System2DCalculations::BerryCurvatureFromWilsonLoop(double kx, double ky, double dk)
{
    double hdk = dk / 2.0;

    Eigen::MatrixXcd Umm = eigenvecs(kx - hdk, ky - hdk);
    Eigen::MatrixXcd Ump = eigenvecs(kx - hdk, ky + hdk);
    Eigen::MatrixXcd Upm = eigenvecs(kx + hdk, ky - hdk);
    Eigen::MatrixXcd Upp = eigenvecs(kx + hdk, ky + hdk);

    Eigen::MatrixXcd W = Eigen::MatrixXcd::Identity(_sys.n_bands, _sys.n_bands);

    _SVD.compute(Umm.leftCols(_sys.n_bands).adjoint() * Upm.leftCols(_sys.n_bands));
    W *= _SVD.matrixU() * _SVD.matrixV().adjoint();

    _SVD.compute(Upm.leftCols(_sys.n_bands).adjoint() * Upp.leftCols(_sys.n_bands));
    W *= _SVD.matrixU() * _SVD.matrixV().adjoint();

    _SVD.compute(Upp.leftCols(_sys.n_bands).adjoint() * Ump.leftCols(_sys.n_bands));
    W *= _SVD.matrixU() * _SVD.matrixV().adjoint();

    _SVD.compute(Ump.leftCols(_sys.n_bands).adjoint() * Umm.leftCols(_sys.n_bands));
    W *= _SVD.matrixU() * _SVD.matrixV().adjoint();

    double BC = std::arg(W.determinant()) / (dk * dk);

    return BC;
}

// bledy numerycznie wywalaja wyniki
double System2DCalculations::MatrixBerryCurvatureTrace(double kx, double ky, double dk)
{
    double hdk = dk / 2.0;

    Eigen::MatrixXcd Uxm = eigenvecs(kx - hdk, ky);
    Eigen::MatrixXcd Uxp = eigenvecs(kx + hdk, ky);
    Eigen::MatrixXcd Uym = eigenvecs(kx, ky - hdk);
    Eigen::MatrixXcd Uyp = eigenvecs(kx, ky + hdk);

    Eigen::MatrixXcd Umm = eigenvecs(kx - hdk, ky - hdk);
    Eigen::MatrixXcd Ump = eigenvecs(kx - hdk, ky + hdk);
    Eigen::MatrixXcd Upm = eigenvecs(kx + hdk, ky - hdk);
    Eigen::MatrixXcd Upp = eigenvecs(kx + hdk, ky + hdk);

    _SAES.compute(_sys.HBdG(kx, ky));

    Eigen::MatrixXcd Ax = 1i * _SAES.eigenvectors().leftCols(_sys.n_bands).adjoint() * ((Uxp - Uxm).leftCols(_sys.n_bands) / dk);
    Eigen::MatrixXcd Ay = 1i * _SAES.eigenvectors().leftCols(_sys.n_bands).adjoint() * ((Uyp - Uym).leftCols(_sys.n_bands) / dk);

    Eigen::MatrixXcd Ay_xm = 1i * Uxm.leftCols(_sys.n_bands).adjoint() * ((Ump - Umm).leftCols(_sys.n_bands) / dk);
    Eigen::MatrixXcd Ay_xp = 1i * Uxp.leftCols(_sys.n_bands).adjoint() * ((Upp - Upm).leftCols(_sys.n_bands) / dk);

    Eigen::MatrixXcd Ax_ym = 1i * Uym.leftCols(_sys.n_bands).adjoint() * ((Upm - Umm).leftCols(_sys.n_bands) / dk);
    Eigen::MatrixXcd Ax_yp = 1i * Uyp.leftCols(_sys.n_bands).adjoint() * ((Upp - Ump).leftCols(_sys.n_bands) / dk);

    Eigen::MatrixXcd B = (Ay_xp - Ay_xm) / dk - (Ax_yp - Ax_ym) / dk + 1i * (Ax * Ay - Ay * Ax); // bez komutatora, bo i tak trace bierzemy

    double BC = B.real().trace();

    // return std::abs(BC) < 1e6 ? BC : 0.0;
    return BC;
}

Eigen::VectorXd System2DCalculations::WilsonLoopSpectrum(std::size_t n, int axis, double k0)
{
    Eigen::VectorXd kx_vec(n);
    Eigen::VectorXd ky_vec(n);

    switch (axis)
    {
    case 0:
        kx_vec = Eigen::VectorXd::Constant(n, k0);
        ky_vec = Eigen::VectorXd::LinSpaced(n, -M_PI, M_PI);
        break;
    case 1:
        kx_vec = Eigen::VectorXd::LinSpaced(n, -M_PI, M_PI);
        ky_vec = Eigen::VectorXd::Constant(n, k0);
        break;
    default:
        throw std::invalid_argument("Invalid axis");
    }

    Eigen::MatrixXcd W = Eigen::MatrixXcd::Identity(_sys.n_bands, _sys.n_bands);

    Eigen::MatrixXcd U0 = eigenvecs(kx_vec(0), ky_vec(0));
    Eigen::MatrixXcd Ui = U0;
    Eigen::MatrixXcd Uip;

    Eigen::MatrixXcd *Ui_ptr = &Ui;
    Eigen::MatrixXcd *Uip_ptr = &Uip;

    for (std::size_t i = 1; i < n - 1; ++i)
    {
        *Uip_ptr = eigenvecs(kx_vec(i), ky_vec(i));
        _SVD.compute((*Ui_ptr).leftCols(_sys.n_bands).adjoint() * (*Uip_ptr).leftCols(_sys.n_bands));
        W *= _SVD.matrixU() * _SVD.matrixV().adjoint();

        std::swap(Ui_ptr, Uip_ptr);
    }

    _SVD.compute(_SAES.eigenvectors().leftCols(_sys.n_bands).adjoint() * U0.leftCols(_sys.n_bands));
    W *= _SVD.matrixU() * _SVD.matrixV().adjoint();

    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> CES(W, Eigen::EigenvaluesOnly);
    Eigen::VectorXd phases = CES.eigenvalues().cwiseArg();

    std::sort(phases.data(), phases.data() + phases.size());

    return phases;
}

Eigen::VectorXd System2DCalculations::ChernNumbersUsingAbelianBerryCurvature(std::size_t n_dense, std::size_t n_sparse, double k_val)
{
    Eigen::VectorXd k_vec = generate_k_vec(n_dense, n_sparse, k_val);

    Eigen::VectorXd phases = Eigen::VectorXd::Zero(2 * _sys.n_bands);

    for (std::size_t i = 0; i < k_vec.size() - 1; ++i)
    {
        for (std::size_t j = 0; j < k_vec.size() - 1; ++j)
        {
            phases += AbelianBerryCurvature((k_vec(i) + k_vec(i + 1)) / 2.0, (k_vec(j) + k_vec(j + 1)) / 2.0) * (k_vec(i + 1) - k_vec(i)) * (k_vec(j + 1) - k_vec(j));
        }
    }

    return phases / (2.0 * M_PI);
}

double System2DCalculations::ChernNumberUsingBerryCurvatureFromWilsonLoop(std::size_t n_dense, std::size_t n_sparse, double k_val)
{
    Eigen::VectorXd k_vec = generate_k_vec(n_dense, n_sparse, k_val);

    double phase = 0.0;

    for (std::size_t i = 0; i < k_vec.size() - 1; ++i)
    {
        for (std::size_t j = 0; j < k_vec.size() - 1; ++j)
        {
            phase += BerryCurvatureFromWilsonLoop((k_vec(i) + k_vec(i + 1)) / 2.0, (k_vec(j) + k_vec(j + 1)) / 2.0) * (k_vec(i + 1) - k_vec(i)) * (k_vec(j + 1) - k_vec(j));
        }
    }

    return phase / (2.0 * M_PI);
}

double System2DCalculations::ChernNumberUsingWilsonLoop(std::size_t n_dense, std::size_t n_sparse, double k_val)
{
    Eigen::VectorXd kx = generate_k_vec(n_dense, n_sparse, k_val);
    Eigen::VectorXd ky = kx;

    std::size_t rows = kx.size();
    std::size_t cols = ky.size();

    Eigen::VectorX<Eigen::MatrixXcd> i_evecs(cols);  // Current row eigenvectors
    Eigen::VectorX<Eigen::MatrixXcd> ni_evecs(cols); // Next row eigenvectors

    Eigen::VectorX<Eigen::MatrixXcd> *i_evecs_ptr = &i_evecs;
    Eigen::VectorX<Eigen::MatrixXcd> *ni_evecs_ptr = &ni_evecs;

    // Compute eigenvectors for the first row
    for (auto j = 0; j < cols; ++j)
    {
        i_evecs(j) = eigenvecs(kx(0), ky(j));
    }

    double CN = 0.0;
    Eigen::MatrixXcd Wcw(_sys.n_bands, _sys.n_bands); // W_{clockwise}

    for (auto i = 0; i < rows - 1; ++i)
    {
        (*ni_evecs_ptr)(0) = eigenvecs(kx(i + 1), ky(0)); // contains eigenvectors and eigenvalues for next i, j=0

        for (auto j = 0; j < cols - 1; ++j)
        {
            (*ni_evecs_ptr)(j + 1) = eigenvecs(kx(i + 1), ky(j + 1)); // contains eigenvectors and eigenvalues for next i, next j

            _SVD.compute((*i_evecs_ptr)(j).leftCols(_sys.n_bands).adjoint() * (*i_evecs_ptr)(j + 1).leftCols(_sys.n_bands));
            Wcw = _SVD.matrixU() * _SVD.matrixV().adjoint();

            _SVD.compute((*i_evecs_ptr)(j + 1).leftCols(_sys.n_bands).adjoint() * (*ni_evecs_ptr)(j + 1).leftCols(_sys.n_bands));
            Wcw *= _SVD.matrixU() * _SVD.matrixV().adjoint();

            _SVD.compute((*ni_evecs_ptr)(j + 1).leftCols(_sys.n_bands).adjoint() * (*ni_evecs_ptr)(j).leftCols(_sys.n_bands));
            Wcw *= _SVD.matrixU() * _SVD.matrixV().adjoint();

            _SVD.compute((*ni_evecs_ptr)(j).leftCols(_sys.n_bands).adjoint() * (*i_evecs_ptr)(j).leftCols(_sys.n_bands));
            Wcw *= _SVD.matrixU() * _SVD.matrixV().adjoint();

            CN += std::arg(Wcw.determinant());
        }

        std::swap(i_evecs_ptr, ni_evecs_ptr);
    }

    return CN / (2.0 * M_PI);
}

double System2DCalculations::ChernNumberUsingWilsonLoop_discrete_kx(std::size_t n_dense, std::size_t n_sparse, double k_val, std::size_t n_kx)
{
    Eigen::VectorXd ky = generate_k_vec(n_dense, n_sparse, k_val);

    std::size_t half_sc_bands = _sys.n_bands * n_kx;

    Eigen::MatrixXcd W = Eigen::MatrixXcd::Identity(half_sc_bands, half_sc_bands);

    Eigen::MatrixXcd U0 = eigenvecs_discrete_kx(ky(0), n_kx);
    Eigen::MatrixXcd Ui = U0;
    Eigen::MatrixXcd Uip;

    Eigen::MatrixXcd *Ui_ptr = &Ui;
    Eigen::MatrixXcd *Uip_ptr = &Uip;

    for (std::size_t i = 1; i < ky.size() - 1; ++i)
    {
        *Uip_ptr = eigenvecs_discrete_kx(ky(i), n_kx);
        _SVD.compute(Ui_ptr->leftCols(half_sc_bands).adjoint() * Uip_ptr->leftCols(half_sc_bands));
        W *= _SVD.matrixU() * _SVD.matrixV().adjoint();

        std::swap(Ui_ptr, Uip_ptr);
    }

    _SVD.compute(_SAES.eigenvectors().leftCols(half_sc_bands).adjoint() * U0.leftCols(half_sc_bands));
    W *= _SVD.matrixU() * _SVD.matrixV().adjoint();

    return std::arg(W.determinant()) / (2.0 * M_PI);
}

double System2DCalculations::ChernNumberUsingWilsonLoop_discrete_ky(std::size_t n_dense, std::size_t n_sparse, double k_val, std::size_t n_ky)
{
    Eigen::VectorXd kx = generate_k_vec(n_dense, n_sparse, k_val);

    std::size_t half_sc_bands = _sys.n_bands * n_ky;

    Eigen::MatrixXcd W = Eigen::MatrixXcd::Identity(half_sc_bands, half_sc_bands);

    Eigen::MatrixXcd U0 = eigenvecs_discrete_ky(kx(0), n_ky);
    Eigen::MatrixXcd Ui = U0;
    Eigen::MatrixXcd Uip;

    Eigen::MatrixXcd *Ui_ptr = &Ui;
    Eigen::MatrixXcd *Uip_ptr = &Uip;

    for (std::size_t i = 1; i < kx.size() - 1; ++i)
    {
        *Uip_ptr = eigenvecs_discrete_ky(kx(i), n_ky);
        _SVD.compute(Ui_ptr->leftCols(half_sc_bands).adjoint() * Uip_ptr->leftCols(half_sc_bands));
        W *= _SVD.matrixU() * _SVD.matrixV().adjoint();

        std::swap(Ui_ptr, Uip_ptr);
    }

    _SVD.compute(_SAES.eigenvectors().leftCols(half_sc_bands).adjoint() * U0.leftCols(half_sc_bands));
    W *= _SVD.matrixU() * _SVD.matrixV().adjoint();

    return std::arg(W.determinant()) / (2.0 * M_PI);
}

std::vector<std::vector<Point2D>> System2DCalculations::FSContours(double E, double dk, double eps, double kx_min, double kx_max, std::size_t n_kx)
{
    std::vector<std::vector<Point2D>> FS_contours;

    // search for kx such that E is an eigenvalue of H(kx, 0)
    Eigen::VectorXd kx_vec = Eigen::VectorXd::LinSpaced(n_kx, kx_min, kx_max);

    _SAES.compute(_sys.Hk(kx_vec(0), 0.0), Eigen::EigenvaluesOnly);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> SAES_next;

    std::vector<std::size_t> kx_idx; // indices of kx_vec where the contour is found
    std::vector<std::size_t> E_idx;  // indices of bands creating the contour

    for (auto ie = 0; ie < _sys.n_bands; ++ie) // contour index
    {
        for (auto i = (kx_idx.empty() ? 0 : kx_idx.back()); i < kx_vec.size() - 1; ++i)
        {
            switch (i % 2)
            {
            case 0:
                SAES_next.compute(_sys.Hk(kx_vec(i + 1), 0.0), Eigen::EigenvaluesOnly);
                break;
            case 1:
                SAES_next.compute(_sys.Hk(kx_vec(i + 1), 0.0), Eigen::EigenvaluesOnly);
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

        double El = _SAES.compute(_sys.Hk(kxl, 0.0), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);
        double Er = _SAES.compute(_sys.Hk(kxr, 0.0), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);

        double kx_exact = 0.0;
        double E_exact = 0.0;

        do
        {
            kx_exact = (El * kxr - Er * kxl) / (El - Er);
            E_exact = _SAES.compute(_sys.Hk(kx_exact, 0.0), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);

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

            El = _SAES.compute(_sys.Hk(kl.x(), kl.y()), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);
            Er = _SAES.compute(_sys.Hk(kr.x(), kr.y()), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);

            // find next k - regula falsi with phi
            do
            {
                phi_next = (El * phir - Er * phil) / (El - Er);
                R_next = R_mat(phi_next);
                k_next = contour.end()[-1] + R_next * k_diff;
                E_next = _SAES.compute(_sys.Hk(k_next.x(), k_next.y()), Eigen::EigenvaluesOnly).eigenvalues()(E_idx[ic]);

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
