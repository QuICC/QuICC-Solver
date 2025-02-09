#pragma once

#include <cuda/std/complex>

namespace QuICC {
namespace Blas {
namespace Cuda {

using namespace QuICC::Memory;

namespace Naive {

/// cuda kernel naive matmul
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

namespace Blocked
{

/// cuda kernel blocked matmul
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


namespace Vad
{

/// cuda kernel vad matmul
/// joint tiling register / shared memory implementation (Volkov and Demmel)
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
        auto kk = threadIdx.x % KK;
        auto nn = threadIdx.x / KK;
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

} // namespace Vad


} // namespace Cuda
} // namespace Blas
} // namespace QuICC

#undef LIBCUDACXX_ENABLE_SIMPLIFIED_COMPLEX_OPERATIONS
