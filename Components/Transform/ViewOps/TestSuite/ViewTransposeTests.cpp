#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "ViewOps/Transpose/Op.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"

using namespace QuICC::Memory;
using namespace QuICC::View;

TEST_CASE("Serial DCCSC3D to DCCSC3D 201", "SerialDCCSC3DtoDCCSC3D201")
{
   // FFT out -> AL in
   // Full data and type
   constexpr size_t M = 4;
   constexpr size_t N = 2;
   constexpr size_t K = 2;

   constexpr size_t S = M * N * K;

   // clang-format off
   std::array<double, S> dataIn = {
      /*k0*/ 1, 2, 3, 4,
      /*k0*/ 5, 6, 7, 8,
      /*k1*/ 9, 10, 11, 12,
      /*k1*/ 13, 14, 15, 16};
   // clang-format on

   // perm = [2 0 1] -> N K M
   std::array<double, S> dataOut;

   // view
   constexpr std::uint32_t rank = 3;
   std::array<std::uint32_t, rank> dimensionsIn{M, N, K};
   std::array<std::uint32_t, rank> dimensionsOut{N, K, M};
   // skip setting up pointers and indices
   // they are not used in the serial aka dense tranpose
   std::array<std::vector<std::uint32_t>, rank> pointers = {{{}, {}, {}}};
   std::array<std::vector<std::uint32_t>, rank> indices = {{{}, {}, {}}};
   using inTy = DCCSC3D;
   using outTy = DCCSC3D;
   View<double, inTy> viewIn(dataIn, dimensionsIn, pointers, indices);
   View<double, outTy> viewOut(dataOut, dimensionsOut, pointers, indices);

   // Transpose op
   using namespace QuICC::Transpose::Cpu;
   using namespace QuICC::Transpose;
   auto transposeOp =
      std::make_unique<Op<View<double, outTy>, View<double, inTy>, p201_t>>();

   transposeOp->apply(viewOut, viewIn);

   // check
   for (std::uint64_t k = 0; k < K; ++k)
   {
      for (std::uint64_t n = 0; n < N; ++n)
      {
         for (std::uint64_t m = 0; m < M; ++m)
         {
            auto mnk = m + n * M + k * M * N;
            auto nkm = n + k * N + m * K * N;
            CHECK(viewIn[mnk] == viewOut[nkm]);
         }
      }
   }
}

TEST_CASE("Serial DCCSC3D to DCCSC3D 120", "SerialDCCSC3DtoDCCSC3D120")
{
   // FFT out -> AL in
   // Full data and type
   constexpr size_t M = 4;
   constexpr size_t N = 2;
   constexpr size_t K = 2;

   constexpr size_t S = M * N * K;
   // clang-format off
   std::array<double, S> dataIn = {
      /*k0*/ 1, 2, 3, 4,
      /*k0*/ 5, 6, 7, 8,
      /*k1*/ 9, 10, 11, 12,
      /*k1*/ 13, 14, 15, 16};
   // clang-format on

   // perm = [1 2 0] -> K M N
   std::array<double, S> dataOut;

   // view
   constexpr std::uint32_t rank = 3;
   std::array<std::uint32_t, rank> dimensionsIn{M, N, K};
   std::array<std::uint32_t, rank> dimensionsOut{K, M, N};
   // skip setting up pointers and indices
   // they are not used in the serial aka dense tranpose
   std::array<std::vector<std::uint32_t>, rank> pointers = {{{}, {}, {}}};
   std::array<std::vector<std::uint32_t>, rank> indices = {{{}, {}, {}}};
   using inTy = DCCSC3D;
   using outTy = DCCSC3D;
   View<double, inTy> viewIn(dataIn, dimensionsIn, pointers, indices);
   View<double, outTy> viewOut(dataOut, dimensionsOut, pointers, indices);

   // Transpose op
   using namespace QuICC::Transpose::Cpu;
   using namespace QuICC::Transpose;
   auto transposeOp =
      std::make_unique<Op<View<double, outTy>, View<double, inTy>, p120_t>>();

   transposeOp->apply(viewOut, viewIn);

   // check
   for (std::uint64_t k = 0; k < K; ++k)
   {
      for (std::uint64_t n = 0; n < N; ++n)
      {
         for (std::uint64_t m = 0; m < M; ++m)
         {
            auto mnk = m + n * M + k * M * N;
            auto kmn = k + m * K + n * K * M;
            CHECK(viewIn[mnk] == viewOut[kmn]);
         }
      }
   }
}

TEST_CASE("Serial S1CLCSC3D to DCCSC3D 201", "SerialS1CLCSC3DtoDCCSC3D201")
{
   // AL out -> JW in
   // Full data and type
   constexpr size_t M = 4;
   constexpr size_t N = 2;
   constexpr size_t K = 3;

   constexpr size_t S = (M + (M - 1) + (M - 2)) * N;
   // clang-format off
   std::array<double, S> dataIn = {
      /*k0*/ 1, 2, 3, 4,
      /*k0*/ 5, 6, 7, 8,
      /*k1*/ 9, 10, 11,
      /*k1*/ 12, 13, 14,
      /*k2*/ 15, 16,
      /*k2*/ 17, 18
   };

   // perm = [2 0 1] -> N K M
   std::array<double, S> dataRef = {
      /*m0*/ 1, 5,
      /*m1*/ 2, 6,
      /*m1*/ 9, 12,
      /*m2*/ 3, 7,
      /*m2*/ 10, 13,
      /*m2*/ 15, 17,
      /*m3*/ 4, 8,
      /*m3*/ 11, 14,
      /*m3*/ 16, 18
   };
   // clang-format on

   std::array<double, S> dataOut;

   // view
   constexpr std::uint32_t rank = 3;
   std::array<std::uint32_t, rank> dimensionsIn{M, N, K};
   std::array<std::uint32_t, rank> dimensionsOut{N, K, M};
   // skip setting up pointers and indices
   // they are not used in the serial aka dense tranpose
   std::array<std::vector<std::uint32_t>, rank> pointers = {{{}, {}, {}}};
   std::array<std::vector<std::uint32_t>, rank> indices = {{{}, {}, {}}};
   using inTy = S1CLCSC3D;
   using outTy = DCCSC3D;
   View<double, inTy> viewIn(dataIn, dimensionsIn, pointers, indices);
   View<double, outTy> viewOut(dataOut, dimensionsOut, pointers, indices);

   // Transpose op
   using namespace QuICC::Transpose::Cpu;
   using namespace QuICC::Transpose;
   auto transposeOp =
      std::make_unique<Op<View<double, outTy>, View<double, inTy>, p201_t>>();

   transposeOp->apply(viewOut, viewIn);

   // check
   for (std::uint64_t s = 0; s < S; ++s)
   {
      CHECK(dataRef[s] == viewOut[s]);
   }
}

TEST_CASE("Serial DCCSC3D to S1CLCSC3D 120", "SerialDCCSC3DtoS1CLCSC3D120")
{
   // JW out -> AL in
   // Full data and type
   constexpr size_t M = 4;
   constexpr size_t N = 2;
   constexpr size_t K = 3;

   constexpr size_t S = (M + (M - 1) + (M - 2)) * N;

   // clang-format off
   // N K M
   std::array<double, S> dataIn = {
      /*m0*/ 1, 5,
      /*m1*/ 2, 6,
      /*m1*/ 9, 12,
      /*m2*/ 3, 7,
      /*m2*/ 10, 13,
      /*m2*/ 15, 17,
      /*m3*/ 4, 8,
      /*m3*/ 11, 14,
      /*m3*/ 16, 18
   };

   // perm = [1 2 0] -> M N K
   std::array<double, S> dataRef = {
      /*k0*/ 1, 2, 3, 4,
      /*k0*/ 5, 6, 7, 8,
      /*k1*/ 9, 10, 11,
      /*k1*/ 12, 13, 14,
      /*k2*/ 15, 16,
      /*k2*/ 17, 18
   };
   // clang-format on

   std::array<double, S> dataOut;

   // view
   constexpr std::uint32_t rank = 3;
   std::array<std::uint32_t, rank> dimensionsIn{N, K, M};
   std::array<std::uint32_t, rank> dimensionsOut{M, N, K};
   // skip setting up pointers and indices
   // they are not used in the serial aka dense tranpose
   std::array<std::vector<std::uint32_t>, rank> pointers = {{{}, {}, {}}};
   std::array<std::vector<std::uint32_t>, rank> indices = {{{}, {}, {}}};
   using inTy = DCCSC3D;
   using outTy = S1CLCSC3D;
   View<double, inTy> viewIn(dataIn, dimensionsIn, pointers, indices);
   View<double, outTy> viewOut(dataOut, dimensionsOut, pointers, indices);

   // Transpose op
   using namespace QuICC::Transpose::Cpu;
   using namespace QuICC::Transpose;
   auto transposeOp =
      std::make_unique<Op<View<double, outTy>, View<double, inTy>, p120_t>>();

   transposeOp->apply(viewOut, viewIn);

   // check
   for (std::uint64_t s = 0; s < S; ++s)
   {
      CHECK(dataRef[s] == viewOut[s]);
   }
}
