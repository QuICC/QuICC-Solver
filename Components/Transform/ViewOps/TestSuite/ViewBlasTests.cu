#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <complex>
#include <cuda/std/complex>

#include "ViewOps/Blas/Cuda/Gemm.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"

/// @brief simple random number generator for floats
/// @tparam T
/// @return a T between -1 and 1
template <class T> inline T randf()
{
   return 2.0 * static_cast<T>(std::rand()) / static_cast<T>(RAND_MAX) - 1.0;
}

enum class memlay
{
   AcmBrmCrm = 0,
};

/// @brief Naive matmul c = a*b to check results
/// @tparam TA
/// @tparam TB
/// @tparam TC
/// @tparam MEM
/// @param c
/// @param a
/// @param b
/// @param M
/// @param K
/// @param N
template <class TC, class TA, class TB, memlay MEM>
void cpu_naive_gemm(TC* c, TA* a, TB* b, const std::size_t M,
   const std::size_t K, const std::size_t N)
{
   for (std::size_t m = 0; m < M; ++m)
   {
      for (std::size_t n = 0; n < N; ++n)
      {
         for (std::size_t k = 0; k < K; ++k)
         {
            if constexpr (MEM == memlay::AcmBrmCrm)
            {
               // a column major
               auto mk = m + k * M;
               // b row major
               auto kn = k * N + n;
               // c row major
               auto mn = m * N + n;
               c[mn] += a[mk] * b[kn];
            }
         }
      }
   }
}


template <class T, unsigned int Nthreads, unsigned int CoarseFactor>
__global__ void wrapperVadRegA(cuda::std::complex<T>* c, T* a,
   cuda::std::complex<T>* b, unsigned int M, unsigned int K, unsigned int N)
{
   QuICC::Blas::Cuda::VadRegA::matmul<T, T, Nthreads, CoarseFactor>(c, a, b, M, K, N, 1.0);
}


TEST_CASE("Mixed GEMM VAD", "MixedGEMMVAD")
{
   constexpr unsigned int M = 1 << 8;
   constexpr unsigned int N = 1 << 8;
   constexpr unsigned int K = 1 << 8;

   std::size_t sizeA = M * K * sizeof(double);
   std::size_t sizeB = K * N * sizeof(std::complex<double>);
   std::size_t sizeC = M * N * sizeof(std::complex<double>);

   std::array<double, M * K> a_h;
   std::array<std::complex<double>, K * N> b_h;
   std::array<std::complex<double>, M * N> c_h;
   std::array<std::complex<double>, M * N> c_r;
   double* a_d;
   cuda::std::complex<double>*b_d, *c_d;

   // Allocate on device
   cudaErrChk(cudaMalloc(&a_d, sizeA));
   cudaErrChk(cudaMalloc(&b_d, sizeB));
   cudaErrChk(cudaMalloc(&c_d, sizeC));

   // Initialize matrices on the host
   for (std::size_t i = 0; i < M * K; ++i)
   {
      a_h[i] = randf<double>();
   }
   for (std::size_t i = 0; i < K * N; ++i)
   {
      b_h[i] = std::complex<double>(randf<double>(), randf<double>());
   }
   for (std::size_t i = 0; i < M * N; ++i)
   {
      c_h[i] = 0.0;
      c_r[i] = 0.0;
   }

   // CPU -> GPU
   cudaErrChk(cudaMemcpy(a_d, a_h.data(), sizeA, cudaMemcpyHostToDevice));
   cudaErrChk(cudaMemcpy(b_d, b_h.data(), sizeB, cudaMemcpyHostToDevice));
   cudaErrChk(cudaMemcpy(c_d, c_h.data(), sizeC, cudaMemcpyHostToDevice));

   // Run kernel on 1M elements on the GPU
   dim3 blockSize;
   dim3 numBlocks;

   constexpr auto layout = memlay::AcmBrmCrm;
   constexpr unsigned int numThreads = 64;
   constexpr unsigned int tileSizeN = 16; // coarsening factor
   blockSize.x = numThreads;
   blockSize.y = 1;
   blockSize.z = 1;
   numBlocks.x = (M + numThreads - 1) / numThreads;
   numBlocks.y = (N + tileSizeN - 1) / tileSizeN;
   numBlocks.z = 1;

   /// vad
   wrapperVadRegA<double, numThreads, tileSizeN><<<numBlocks, blockSize>>>(c_d, a_d, b_d, M, K, N);

   // GPU -> CPU
   cudaErrChk(cudaMemcpy(c_h.data(), c_d, sizeC, cudaMemcpyDeviceToHost));

   /// check
   cpu_naive_gemm<std::complex<double>, double, std::complex<double>, layout>( c_r.data(), a_h.data(), b_h.data(), M, K, N);

   double eps = 1.e-10;
   for (std::uint64_t i = 0; i < c_r.size(); ++i)
   {
      CHECK(std::abs(c_r[i].real() - c_h[i].real()) < eps );
      CHECK(std::abs(c_r[i].imag() - c_h[i].imag()) < eps );
   }

   // dealloc device
   cudaErrChk(cudaFree(a_d));
   cudaErrChk(cudaFree(b_d));
   cudaErrChk(cudaFree(c_d));
}

template <class T, unsigned int Nthreads, unsigned int CoarseFactor>
__global__ void wrapperVadRegB(cuda::std::complex<T>* c, T* a,
   cuda::std::complex<T>* b, unsigned int M, unsigned int K, unsigned int N)
{
   QuICC::Blas::Cuda::VadRegB::matmul<T, T, Nthreads, CoarseFactor>(c, a, b, M, K, N, 1.0);
}


TEST_CASE("Mixed GEMM VAD reg b", "MixedGEMMVADRegB")
{
   constexpr unsigned int M = 1 << 8;
   constexpr unsigned int N = 1 << 8;
   constexpr unsigned int K = 1 << 8;

   std::size_t sizeA = M * K * sizeof(double);
   std::size_t sizeB = K * N * sizeof(std::complex<double>);
   std::size_t sizeC = M * N * sizeof(std::complex<double>);

   std::array<double, M * K> a_h;
   std::array<std::complex<double>, K * N> b_h;
   std::array<std::complex<double>, M * N> c_h;
   std::array<std::complex<double>, M * N> c_r;
   double* a_d;
   cuda::std::complex<double>*b_d, *c_d;

   // Allocate on device
   cudaErrChk(cudaMalloc(&a_d, sizeA));
   cudaErrChk(cudaMalloc(&b_d, sizeB));
   cudaErrChk(cudaMalloc(&c_d, sizeC));

   // Initialize matrices on the host
   for (std::size_t i = 0; i < M * K; ++i)
   {
      a_h[i] = randf<double>();
   }
   for (std::size_t i = 0; i < K * N; ++i)
   {
      b_h[i] = std::complex<double>(randf<double>(), randf<double>());
   }
   for (std::size_t i = 0; i < M * N; ++i)
   {
      c_h[i] = 0.0;
      c_r[i] = 0.0;
   }

   // CPU -> GPU
   cudaErrChk(cudaMemcpy(a_d, a_h.data(), sizeA, cudaMemcpyHostToDevice));
   cudaErrChk(cudaMemcpy(b_d, b_h.data(), sizeB, cudaMemcpyHostToDevice));
   cudaErrChk(cudaMemcpy(c_d, c_h.data(), sizeC, cudaMemcpyHostToDevice));

   // Run kernel on 1M elements on the GPU
   dim3 blockSize;
   dim3 numBlocks;

   constexpr auto layout = memlay::AcmBrmCrm;
   constexpr unsigned int numThreads = 64;
   constexpr unsigned int tileSizeM = 16; // coarsening factor
   blockSize.x = numThreads;
   blockSize.y = 1;
   blockSize.z = 1;
   numBlocks.x = (M + tileSizeM - 1) / tileSizeM;
   numBlocks.y = (N + numThreads - 1) / numThreads;
   numBlocks.z = 1;

   /// vad
   wrapperVadRegB<double, numThreads, tileSizeM><<<numBlocks, blockSize>>>(c_d, a_d, b_d, M, K, N);

   // GPU -> CPU
   cudaErrChk(cudaMemcpy(c_h.data(), c_d, sizeC, cudaMemcpyDeviceToHost));

   /// check
   cpu_naive_gemm<std::complex<double>, double, std::complex<double>, layout>( c_r.data(), a_h.data(), b_h.data(), M, K, N);

   double eps = 1.e-10;
   for (std::uint64_t i = 0; i < c_r.size(); ++i)
   {
      CHECK(std::abs(c_r[i].real() - c_h[i].real()) < eps );
      CHECK(std::abs(c_r[i].imag() - c_h[i].imag()) < eps );
   }

   // dealloc device
   cudaErrChk(cudaFree(a_d));
   cudaErrChk(cudaFree(b_d));
   cudaErrChk(cudaFree(c_d));
}
