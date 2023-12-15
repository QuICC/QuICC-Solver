/**
 * @file Gemm.hpp
 * @brief Different implementations of Cpu Gemm operations on Views
 */
#pragma once

// External includes
//
#include <Eigen/Core>

// Project includes
//
#include "View/View.hpp"

namespace QuICC {
/// @brief namespace for Blas type operations
namespace Blas {
/// @brief namespace for cpu backends
namespace Cpu {

/// @brief namespace for Naive implementations
namespace Naive {

/// @brief Naive matmul with mixed types
/// @tparam TA row major
/// @tparam TB column major
/// @tparam TC column major
/// @tparam Talpha
/// @param C
/// @param A
/// @param B
/// @param alpha
template <class TA, class TB, class TC, class Talpha>
inline void matmul(View::View<TC, View::dense2D>& C,
    const View::View<TA, View::dense2DRM>& A,
    const View::View<TB, View::dense2D>& B,
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

/// @brief Naive matmul with mixed types
/// @tparam TA column major
/// @tparam TB row major
/// @tparam TC row major
/// @tparam Talpha
/// @param C
/// @param A
/// @param B
/// @param alpha
template <class TA, class TB, class TC, class Talpha>
inline void matmul(View::View<TC, View::dense2DRM>& C,
    const View::View<TA, View::dense2D>& A,
    const View::View<TB, View::dense2DRM>& B,
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

/// @brief namespace for Eigen backend
namespace Eigen {

/// @brief Eigen matmul with mixed types
/// @tparam TA row major
/// @tparam TB column major
/// @tparam TC column major
/// @tparam Talpha
/// @param C
/// @param A
/// @param B
/// @param alpha
template <class TA, class TB, class TC, class Talpha>
inline void matmul(View::View<TC, View::dense2D>& C,
    const View::View<TA, View::dense2DRM>& A,
    const View::View<TB, View::dense2D>& B,
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

    // Order of operations is important with mixed (complex/real) types!
    eC = eA * eB * alpha;
}

/// @brief Eigen matmul with mixed types
/// @tparam TA column major
/// @tparam TB row major
/// @tparam TC row major
/// @tparam Talpha
/// @param C
/// @param A
/// @param B
/// @param alpha
template <class TA, class TB, class TC, class Talpha>
inline void matmul(View::View<TC, View::dense2DRM>& C,
    const View::View<TA, View::dense2D>& A,
    const View::View<TB, View::dense2DRM>& B,
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

    // Order of operations is important with mixed (complex/real) types!
    eC = eA * eB * alpha;
}

} // namespace Eigen

} // namespace Cpu
} // namespace Blas
} // namespace QuICC
