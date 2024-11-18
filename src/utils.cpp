#include "utils.h"

Eigen::MatrixXcd kron(const Eigen::MatrixXcd &A, const Eigen::MatrixXcd &B)
{
    Eigen::MatrixXcd res(A.rows() * B.rows(), A.cols() * B.cols());

    for (auto i = 0; i < A.rows(); ++i)
    {
        for (auto j = 0; j < A.cols(); ++j)
        {
            res.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i, j) * B;
        }
    }

    return res;
}