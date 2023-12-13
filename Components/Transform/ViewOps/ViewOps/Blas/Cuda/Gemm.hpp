/**
 * @file Gemm.hpp
 * @brief Different implementations of Cuda Gemm operations
 */
#pragma once

// External includes
//
#include <cuda/std/complex>

// Project includes
//

namespace QuICC {
/// @brief namespace for Blas type operations
namespace Blas {
/// @brief namespace for Cuda backends
namespace Cuda {

/// @brief namespace for naive backend
namespace Naive {


/// @brief Cuda kernel for naive matmul with mixed types
/// @tparam T
/// @tparam Talpha
/// @param c complex, row major
/// @param a real, col major
/// @param b complex, row major
/// @param M
/// @param K
/// @param N
/// @param alpha
template <class T, class Talpha>
inline __device__ void matmul(cuda::std::complex<T>* c, const T* a, const cuda::std::complex<T>* b,
    const std::size_t M, const std::size_t K, const std::size_t N,
    const Talpha alpha)
{
    const std::size_t m = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t n = blockIdx.y * blockDim.y + threadIdx.y;

    if(m < M && n < N) {
        cuda::std::complex<T> acc = 0.0;
        auto mn = m*N + n;
        for (std::size_t k = 0; k < K; ++k)
        {
            // row major
            auto kn = k*N + n;
            // col major
            auto mk = m + k*M;
            acc += a[mk] * b[kn];
        }
        c[mn] = alpha * acc;
    }
}

} // namespace Naive

/// @brief namespace for one level blocked backend
namespace Blocked
{

/// @brief Cuda kernel for one level blocked matmul with mixed types
/// @tparam T
/// @tparam Talpha
/// @param c complex, row major
/// @param a real, col major
/// @param b complex, row major
/// @param M
/// @param K
/// @param N
/// @param alpha
template <class T, class Talpha, unsigned int TileSize>
inline __device__ void matmul(cuda::std::complex<T>* c, const T* a, const cuda::std::complex<T>* b,
    const std::size_t M, const std::size_t K, const std::size_t N,
    const Talpha alpha)
{
    // tile size
    constexpr unsigned int tk = TileSize;
    constexpr unsigned int tm = tk;
    constexpr unsigned int tn = tk;

    assert(blockDim.x == tm);
    assert(blockDim.y == tn);
    assert(blockDim.z == 1);

    // local memory
    __shared__ T a_sm[tk][tk];
    __shared__ cuda::std::complex<T> b_sm[tk][tk];

    const auto mm = blockIdx.x*tm;
    const auto nn = blockIdx.y*tn;

    const auto m = threadIdx.y;
    const auto n = threadIdx.x;


    // thread-private accumulator
    cuda::std::complex<T> cp{0.0};
    for (std::size_t kk = 0; kk < K; kk+=tk)
    {
        // copy to local memory (SM)
        // column major
        if (mm + m < M && kk + n < K) {
            a_sm[m][n] = a[mm + m +(kk + n)*M];
        }
        else {
            a_sm[m][n] = 0.0;
        }
        // row major
        if (kk + m < K && nn + n < N) {
            b_sm[m][n] = b[(kk + m)*N + nn + n];
        }
        else {
            b_sm[m][n] = 0.0;
        }
        __syncthreads();

        for (std::size_t k = 0; k < tk; ++k)
        {
            // row major
            cp += a_sm[m][k] * b_sm[k][n];
        }

        __syncthreads();
    }
    if (mm + m < M && nn + n < N) {
        // row major
        c[(mm + m)*N + nn + n] += alpha * cp;
    }
}

} // namespace Blocked

/// @brief namespace for joint tiling and register backend
namespace VadRegA
{

/// @brief Cuda kernel vad matmul
/// joint tiling register / shared memory implementation (Volkov and Demmel)
/// tile a in registers
/// @tparam T
/// @tparam Talpha
/// @param c complex, row major
/// @param a real, col major
/// @param b complex, row major
/// @param M
/// @param K
/// @param N
/// @param alpha
template <class T, class Talpha, unsigned int Nthreads, unsigned int CoarseFactor>
inline __device__ void matmul(cuda::std::complex<T>* c, const T* a, const cuda::std::complex<T>* b,
    const std::size_t M, const std::size_t K, const std::size_t N,
    const Talpha alpha)
{
    // tile size
    constexpr unsigned int ntrd = Nthreads;       // number of threads, min one warp ~64
    constexpr unsigned int MM = ntrd;
    constexpr unsigned int NN = CoarseFactor;     // coarsening factor,
    // number of reuse of reg values, adjust register pressure ~16
    constexpr unsigned int KK = MM / NN;  // steps of outer product, chosen so that tile b_sm is loaded by all threads

    assert(blockDim.x == MM);
    assert(blockDim.y == 1);
    assert(blockDim.z == 1);

    // local memory
    __shared__ cuda::std::complex<T> b_sm[KK*NN];

    // registers
    T a_reg[KK];
    cuda::std::complex<T> c_reg[NN];

    const auto m = blockIdx.x*MM;
    const auto n = blockIdx.y*NN;

    // set tile c to zero
    #pragma unroll
    for (std::uint32_t nn = 0; nn < NN; ++nn) {
        c_reg[nn] = cuda::std::complex<T>{0.0};
    }

    auto mm = threadIdx.x;

    for (std::uint32_t k = 0; k < K; k+=KK)
    {
        // copy b to local memory (SM)
        auto kk = threadIdx.x / NN;
        auto nn = threadIdx.x % NN;
        // row major
        auto kknn = kk*NN + nn;
        if (k + kk < K && n + nn < N) {
            b_sm[kknn] = b[(k + kk)*N + n + nn]; // row major
        }
        else {
            b_sm[kknn] = 0.0;
        }

        // copy a to registers
        #pragma unroll
        for (std::uint32_t kk = 0; kk < KK; ++kk) {
            // col major
            auto mmkk = m + mm + (k + kk)*M;
            if (k + kk < K && m + mm < M) {
                a_reg[kk] = a[mmkk];
            }
            else {
                a_reg[kk] = 0.0;
            }
        }

        __syncthreads();

        // tile c
        #pragma unroll
        for (std::uint32_t nn = 0; nn < NN; ++nn) {
            #pragma unroll
            // KK steps of tile c
            for (std::uint32_t kk = 0; kk < KK; ++kk) {
                // row major
                auto kknn = kk*NN + nn;
                c_reg[nn] += a_reg[kk]*b_sm[kknn];
            }
        }

        __syncthreads();
    }

    // copy back c_reg
    #pragma unroll
    for (std::uint32_t nn = 0; nn < NN; ++nn) {
        if (m + mm < M && n + nn < N) {
            // row major
            c[(mm + m)*N + nn + n] = alpha * c_reg[nn];
        }
    }
}

} // namespace VadRegA

/// @brief namespace for joint tiling and register backend
namespace VadRegB
{

/// @brief Cuda kernel vad matmul
/// joint tiling register / shared memory implementation (Volkov and Demmel)
/// tile b in registers
/// @tparam T
/// @tparam Talpha
/// @param c complex, row major
/// @param a real, col major
/// @param b complex, row major
/// @param M
/// @param K
/// @param N
/// @param alpha
template <class T, class Talpha, unsigned int Nthreads, unsigned int CoarseFactor>
inline __device__ void matmul(cuda::std::complex<T>* c, const T* a, const cuda::std::complex<T>* b,
    const std::size_t M, const std::size_t K, const std::size_t N,
    const Talpha alpha)
{
    // tile size
    constexpr unsigned int ntrd = Nthreads;       // number of threads, min one warp ~64
    constexpr unsigned int MM = CoarseFactor;     // coarsening factor,
    constexpr unsigned int NN = ntrd;
    // number of reuse of reg values, adjust register pressure ~16
    constexpr unsigned int KK = NN / MM;  // steps of outer product, chosen so that tile b_sm is loaded by all threads

    assert(blockDim.x == NN);
    assert(blockDim.y == 1);
    assert(blockDim.z == 1);

    // local memory
    __shared__ T a_sm[MM*KK];

    // registers
    cuda::std::complex<T> b_reg[KK];
    cuda::std::complex<T> c_reg[MM];

    const auto m = blockIdx.x*MM;
    const auto n = blockIdx.y*NN;

    // set tile c to zero
    #pragma unroll
    for (std::uint32_t mm = 0; mm < MM; ++mm) {
        c_reg[mm] = cuda::std::complex<T>{0.0};
    }

    auto nn = threadIdx.x;

    for (std::uint32_t k = 0; k < K; k+=KK)
    {
        // copy a to local memory (SM)
        auto mm = threadIdx.x % MM;
        auto kk = threadIdx.x / MM;
        // col major
        auto mmkk = mm + kk*MM;
        if (k + kk < K && m + mm < M) {
            a_sm[mmkk] = a[m + mm + (k + kk)*M]; // col major
        }
        else {
            a_sm[mmkk] = 0.0;
        }

        // copy b to registers
        #pragma unroll
        for (std::uint32_t kk = 0; kk < KK; ++kk) {
            // row major
            auto kknn = (k + kk)*N + n + nn;
            if (k + kk < K && n + nn < N) {
                b_reg[kk] = b[kknn];
            }
            else {
                b_reg[kk] = 0.0;
            }
        }

        __syncthreads();

        // tile c
        #pragma unroll
        for (std::uint32_t mm = 0; mm < MM; ++mm) {
            #pragma unroll
            // KK steps of tile c
            for (std::uint32_t kk = 0; kk < KK; ++kk) {
                // col major
                auto mmkk = mm + kk*MM;
                c_reg[mm] += a_sm[mmkk]*b_reg[kk];
            }
        }

        __syncthreads();
    }

    // copy back c_reg
    #pragma unroll
    for (std::uint32_t mm = 0; mm < MM; ++mm) {
        if (m + mm < M && n + nn < N) {
            // row major
            c[(mm + m)*N + nn + n] = alpha * c_reg[mm];
        }
    }
}

} // namespace VadRegB


} // namespace Cuda
} // namespace Blas
} // namespace QuICC
