#pragma once

#include <Eigen/Core>

#include "View/View.hpp"

namespace QuICC {
namespace Blas {
namespace Cpu {

using namespace QuICC::Memory;

namespace Naive {

template <class TA, class TB, class TC, class Talpha>
inline void matmul(View<TC, dense2D>& C,
    const View<TA, dense2DRM>& A,
    const View<TB, dense2D>& B,
    const Talpha alpha)
{
    assert(C.dims()[0] == A.dims()[0]);
    assert(C.dims()[1] == B.dims()[1]);
    assert(A.dims()[1] == B.dims()[0]);

    const auto M = C.dims()[0];
    const auto N = C.dims()[1];
    const auto K = A.dims()[1];

    for (std::size_t n = 0; n < N; ++n)
    {
        for (std::size_t m = 0; m < M; ++m)
        {
            TC acc{};
            for (std::size_t k = 0; k < K; ++k)
            {
                // C-B are column major
                // A is row major
                acc += A.data()[k + m*K] * B.data()[k + n*K];
            }
            C.data()[m +n*M] = alpha * acc;
        }
    }
}

template <class TA, class TB, class TC, class Talpha>
inline void matmul(View<TC, dense2DRM>& C,
    const View<TA, dense2D>& A,
    const View<TB, dense2DRM>& B,
    const Talpha alpha)
{
    assert(C.dims()[0] == A.dims()[0]);
    assert(C.dims()[1] == B.dims()[1]);
    assert(A.dims()[1] == B.dims()[0]);

    const auto M = C.dims()[0];
    const auto N = C.dims()[1];
    const auto K = A.dims()[1];

    for (std::size_t n = 0; n < N; ++n)
    {
        for (std::size_t m = 0; m < M; ++m)
        {
            TC acc{};
            for (std::size_t k = 0; k < K; ++k)
            {
                // C-B are row major
                // A is column major
                acc += A.data()[k*M + m] * B.data()[k*N + n];
            }
            C.data()[m*N + n] = alpha * acc;
        }
    }
}

} // namespace Naive

namespace Eigen {

template <class TA, class TB, class TC, class Talpha>
inline void matmul(View<TC, dense2D>& C,
    const View<TA, dense2DRM>& A,
    const View<TB, dense2D>& B,
    const Talpha alpha)
{
    assert(C.dims()[0] == A.dims()[0]);
    assert(C.dims()[1] == B.dims()[1]);
    assert(A.dims()[1] == B.dims()[0]);

    const auto M = C.dims()[0];
    const auto N = C.dims()[1];
    const auto K = A.dims()[1];

    using AMatrixRM = ::Eigen::Matrix<TA, ::Eigen::Dynamic, ::Eigen::Dynamic, ::Eigen::RowMajor>;
    using BMatrixZ = ::Eigen::Matrix<TB, ::Eigen::Dynamic, ::Eigen::Dynamic>;
    using CMatrixZ = ::Eigen::Matrix<TC, ::Eigen::Dynamic, ::Eigen::Dynamic>;

    ::Eigen::Map<AMatrixRM> eA(A.data(), M, K);
    ::Eigen::Map<BMatrixZ> eB(B.data(), K, N);
    ::Eigen::Map<CMatrixZ> eC(C.data(), M, N);

    eC = eA * eB * alpha;
}

template <class TA, class TB, class TC, class Talpha>
inline void matmul(View<TC, dense2DRM>& C,
    const View<TA, dense2D>& A,
    const View<TB, dense2DRM>& B,
    const Talpha alpha)
{
    assert(C.dims()[0] == A.dims()[0]);
    assert(C.dims()[1] == B.dims()[1]);
    assert(A.dims()[1] == B.dims()[0]);

    const auto M = C.dims()[0];
    const auto N = C.dims()[1];
    const auto K = A.dims()[1];

    using AMatrix = ::Eigen::Matrix<TA, ::Eigen::Dynamic, ::Eigen::Dynamic>;
    using BMatrixZRM = ::Eigen::Matrix<TB, ::Eigen::Dynamic, ::Eigen::Dynamic, ::Eigen::RowMajor>;
    using CMatrixZRM = ::Eigen::Matrix<TC, ::Eigen::Dynamic, ::Eigen::Dynamic, ::Eigen::RowMajor>;

    ::Eigen::Map<AMatrix> eA(A.data(), M, K);
    ::Eigen::Map<BMatrixZRM> eB(B.data(), K, N);
    ::Eigen::Map<CMatrixZRM> eC(C.data(), M, N);

    eC = eA * eB * alpha;
}

} // namespace Eigen

} // namespace Cpu
} // namespace Blas
} // namespace QuICC
