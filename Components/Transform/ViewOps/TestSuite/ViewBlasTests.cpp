#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <complex>

#include "View/ViewDense.hpp"
#include "ViewOps/Blas/Cpu/Gemm.hpp"

/// @brief simple random number generator for floats
/// @tparam T
/// @return a T between -1 and 1
template <class T> inline T randf()
{
   return 2.0 * static_cast<T>(std::rand()) / static_cast<T>(RAND_MAX) - 1.0;
}

enum class memlay
{
   ArmBcmCcm = 0,
   AcmBrmCrm = 1,
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
            else if constexpr (MEM == memlay::ArmBcmCcm)
            {
               // a column major
               auto mk = m*K + k;
               // b row major
               auto kn = k + n * K;
               // c row major
               auto mn = m + n * M;
               c[mn] += a[mk] * b[kn];
            }
         }
      }
   }
}

TEST_CASE("Mixed GEMM using Naive ArmBcmCcm", "[MixedGEMMNaiveArmBcmCcm]")
{
   constexpr unsigned int M = 1 << 8;
   constexpr unsigned int N = 1 << 8;
   constexpr unsigned int K = 1 << 8;

   std::array<double, M * K> a_h;
   std::array<std::complex<double>, K * N> b_h;
   std::array<std::complex<double>, M * N> c_h;
   std::array<std::complex<double>, M * N> c_r;

   // Initialize matrices
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

   std::array<std::uint32_t, 2> dimA {M, K};
   std::array<std::uint32_t, 2> dimB {K, N};
   std::array<std::uint32_t, 2> dimC {M, N};

   QuICC::View::View<double, QuICC::View::dense2DRM> vA(a_h, dimA);
   QuICC::View::View<std::complex<double>, QuICC::View::dense2D> vB(b_h, dimB);
   QuICC::View::View<std::complex<double>, QuICC::View::dense2D> vC(c_h, dimC);

   constexpr auto layout = memlay::ArmBcmCcm;

   /// Naive
   QuICC::Blas::Cpu::Naive::matmul(vC, vA, vB, 1.0);

   /// check
   cpu_naive_gemm<std::complex<double>, double, std::complex<double>, layout>( c_r.data(), a_h.data(), b_h.data(), M, K, N);

   double eps = 1.e-10;
   for (std::uint64_t i = 0; i < c_r.size(); ++i)
   {
      CHECK(std::abs(c_r[i].real() - c_h[i].real()) < eps );
      CHECK(std::abs(c_r[i].imag() - c_h[i].imag()) < eps );
   }
}

TEST_CASE("Mixed GEMM using Naive AcmBrmCrm", "[MixedGEMMNaiveAcmBrmCrm]")
{
   constexpr unsigned int M = 1 << 8;
   constexpr unsigned int N = 1 << 8;
   constexpr unsigned int K = 1 << 8;

   std::array<double, M * K> a_h;
   std::array<std::complex<double>, K * N> b_h;
   std::array<std::complex<double>, M * N> c_h;
   std::array<std::complex<double>, M * N> c_r;

   // Initialize matrices
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

   std::array<std::uint32_t, 2> dimA {M, K};
   std::array<std::uint32_t, 2> dimB {K, N};
   std::array<std::uint32_t, 2> dimC {M, N};

   QuICC::View::View<double, QuICC::View::dense2D> vA(a_h, dimA);
   QuICC::View::View<std::complex<double>, QuICC::View::dense2DRM> vB(b_h, dimB);
   QuICC::View::View<std::complex<double>, QuICC::View::dense2DRM> vC(c_h, dimC);

   constexpr auto layout = memlay::AcmBrmCrm;

   /// Naive
   QuICC::Blas::Cpu::Naive::matmul(vC, vA, vB, 1.0);

   /// check
   cpu_naive_gemm<std::complex<double>, double, std::complex<double>, layout>( c_r.data(), a_h.data(), b_h.data(), M, K, N);

   double eps = 1.e-10;
   for (std::uint64_t i = 0; i < c_r.size(); ++i)
   {
      CHECK(std::abs(c_r[i].real() - c_h[i].real()) < eps );
      CHECK(std::abs(c_r[i].imag() - c_h[i].imag()) < eps );
   }
}

TEST_CASE("Mixed GEMM using Eigen AcmBrmCrm", "[MixedGEMMEigenAcmBrmCrm]")
{
   constexpr unsigned int M = 1 << 8;
   constexpr unsigned int N = 1 << 8;
   constexpr unsigned int K = 1 << 8;

   std::array<double, M * K> a_h;
   std::array<std::complex<double>, K * N> b_h;
   std::array<std::complex<double>, M * N> c_h;
   std::array<std::complex<double>, M * N> c_r;

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

   std::array<std::uint32_t, 2> dimA {M, K};
   std::array<std::uint32_t, 2> dimB {K, N};
   std::array<std::uint32_t, 2> dimC {M, N};

   QuICC::View::View<double, QuICC::View::dense2D> vA(a_h, dimA);
   QuICC::View::View<std::complex<double>, QuICC::View::dense2DRM> vB(b_h, dimB);
   QuICC::View::View<std::complex<double>, QuICC::View::dense2DRM> vC(c_h, dimC);

   constexpr auto layout = memlay::AcmBrmCrm;

   /// Eigen
   QuICC::Blas::Cpu::Eigen::matmul(vC, vA, vB, 1.0);

   /// check
   cpu_naive_gemm<std::complex<double>, double, std::complex<double>, layout>( c_r.data(), a_h.data(), b_h.data(), M, K, N);

   double eps = 1.e-10;
   for (std::uint64_t i = 0; i < c_r.size(); ++i)
   {
      CHECK(std::abs(c_r[i].real() - c_h[i].real()) < eps );
      CHECK(std::abs(c_r[i].imag() - c_h[i].imag()) < eps );
   }
}

TEST_CASE("Mixed GEMM using Eigen ArmBcmCcm", "[MixedGEMMEigenArmBcmCcm]")
{
   constexpr unsigned int M = 1 << 8;
   constexpr unsigned int N = 1 << 8;
   constexpr unsigned int K = 1 << 8;

   std::array<double, M * K> a_h;
   std::array<std::complex<double>, K * N> b_h;
   std::array<std::complex<double>, M * N> c_h;
   std::array<std::complex<double>, M * N> c_r;

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

   std::array<std::uint32_t, 2> dimA {M, K};
   std::array<std::uint32_t, 2> dimB {K, N};
   std::array<std::uint32_t, 2> dimC {M, N};

   QuICC::View::View<double, QuICC::View::dense2DRM> vA(a_h, dimA);
   QuICC::View::View<std::complex<double>, QuICC::View::dense2D> vB(b_h, dimB);
   QuICC::View::View<std::complex<double>, QuICC::View::dense2D> vC(c_h, dimC);

   constexpr auto layout = memlay::ArmBcmCcm;

   /// Eigen
   QuICC::Blas::Cpu::Eigen::matmul(vC, vA, vB, 1.0);

   /// check
   cpu_naive_gemm<std::complex<double>, double, std::complex<double>, layout>( c_r.data(), a_h.data(), b_h.data(), M, K, N);

   double eps = 1.e-10;
   for (std::uint64_t i = 0; i < c_r.size(); ++i)
   {
      CHECK(std::abs(c_r[i].real() - c_h[i].real()) < eps );
      CHECK(std::abs(c_r[i].imag() - c_h[i].imag()) < eps );
   }
}
