#define CATCH_CONFIG_RUNNER

#include <catch2/catch.hpp>
#include <chrono>
#include <climits>
#include <regex>
#include <thread>

#include "ViewOps/Transpose/Mpi/Comm.hpp"
extern "C" {
#include <unistd.h>
}

#include "ViewOps/Transpose/Transpose.hpp"
#include "ViewOps/ViewIndexUtils.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#include "TestSuite/ViewMeta.hpp"
#include "Environment/QuICCEnv.hpp"


using namespace QuICC::Transpose::Mpi;


int main(int argc, char** argv)
{
   QuICC::QuICCEnv();

   Catch::Session session; // There must be exactly one instance

   auto returnCode = session.run();

   return returnCode;
}

using namespace QuICC::Memory;
using namespace QuICC::View;
using namespace QuICC::TestSuite;

TEST_CASE("Mpi DCCSC3D to DCCSC3D 210", "MpiDCCSC3DtoDCCSC3D210")
{
   int rank, ranks;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &ranks);

   // FFT out -> AL in
   // Full data and type
   constexpr size_t M = 4;
   constexpr size_t N = 2;
   constexpr size_t K = 2;

   // view
   constexpr std::uint32_t vRank = 3;
   std::array<std::uint32_t, vRank> dimensionsIn{M, N, K};
   std::array<std::uint32_t, vRank> dimensionsOut{N, K, M}; // 2 0 1

   std::vector<double> dataIn;
   std::vector<double> dataOut;
   std::vector<double> dataOutRef;
   ptrAndIdx metaIn;
   ptrAndIdx metaOut;

   if (rank == 0 && ranks == 1)
   {
      constexpr size_t S = M * N * K;
      dataIn = {/*k0*/ 1, 2, 3, 4,
         /*k0*/ 5, 6, 7, 8,
         /*k1*/ 9, 10, 11, 12,
         /*k1*/ 13, 14, 15, 16};

      // perm = [2 0 1] -> N K M
      dataOut.resize(S);

      // Populate meta for fully populated tensor
      // Physical space (Stage::PPP and Stage::MPP)
      metaIn = Index::densePtrAndIdx<DCCSC3D>(dimensionsIn);
      // Populate meta for fully populated tensor
      // AL space (Stage::PPM and Stage::MPM)
      metaOut = Index::densePtrAndIdx<DCCSC3D>(dimensionsOut);
   }
   else if (rank == 0 && ranks == 2)
   {
      metaIn.ptr = {0, 2, 2};
      metaIn.idx = {0, 1};
      dataIn.resize(metaIn.idx.size()*M);
      // perm = [2 0 1] -> N K M
      metaOut.ptr = {0, 2, 4, 4, 4};
      metaOut.idx = {0, 1, 0, 1};
      dataOut.resize(metaOut.idx.size()*N);
      dataOutRef.resize(metaOut.idx.size()*N);
   }
   else if (rank == 1 && ranks == 2)
   {
      metaIn.ptr = {0, 0, 2};
      metaIn.idx = {0, 1};
      dataIn.resize(metaIn.idx.size()*M);
      // perm = [2 0 1] -> N K M
      metaOut.ptr = {0, 0, 0, 2, 4};
      metaOut.idx = {0, 1, 0, 1};
      dataOut.resize(metaOut.idx.size()*N);
      dataOutRef.resize(metaOut.idx.size()*N);
   }
   else if (ranks == 4)
   {
      std::string path = "/home/gcastigl/codes/QuICC/build/Components/Framework/TestSuite/_refdata/Framework/LoadSplitter/WLFl/";
      std::string dist = "Tubular";
      std::string id = "103";
      auto setup = readDimsAndMeta(path, dist, id);

      dimensionsIn = {setup.modsDims[2], setup.physDims[1], setup.physDims[0]};
      dimensionsOut = {setup.physDims[1], setup.physDims[0], setup.physDims[2]};

      metaIn.ptr = setup.metaFT.ptr;
      metaIn.idx = setup.metaFT.idx;
      dataIn.resize(metaIn.idx.size()*dimensionsIn[0]);
      // perm = [2 0 1] MLN -> LNM
      metaOut.ptr = setup.metaAL.ptr;
      metaOut.idx = setup.metaAL.idx;
      dataOut.resize(metaOut.idx.size()*dimensionsOut[0]);
      dataOutRef.resize(metaOut.idx.size()*dimensionsOut[0]);
   }

   // Set meta
   std::array<std::vector<std::uint32_t>, vRank> pointersIn = {
      {{}, metaIn.ptr, {}}};
   std::array<std::vector<std::uint32_t>, vRank> indicesIn = {
      {{}, metaIn.idx, {}}};
   std::array<std::vector<std::uint32_t>, vRank> pointersOut = {
      {{}, metaOut.ptr, {}}};
   std::array<std::vector<std::uint32_t>, vRank> indicesOut = {
      {{}, metaOut.idx, {}}};

   using inTy = DCCSC3D;
   using outTy = DCCSC3D;
   View<double, inTy> viewIn(dataIn, dimensionsIn, pointersIn, indicesIn);
   View<double, outTy> viewOut(dataOut, dimensionsOut, pointersOut, indicesOut);

   // Setup ref data and input data
   using namespace QuICC::Transpose::Mpi;
   using namespace QuICC::Transpose;
   if (ranks > 1)
   {
      auto cooOld = ::QuICC::View::getCoo<View<double, inTy>, p012_t>(viewIn);
      double shift = 2048;
      for (std::size_t i = 0; i < cooOld.size(); ++i)
      {
         dataIn[i] = cooOld[i][0]+cooOld[i][1]*shift+cooOld[i][2]/shift;
      }
      auto cooNew = ::QuICC::View::getCoo<View<double, outTy>, p201_t>(viewOut);
      for (std::size_t i = 0; i < cooNew.size(); ++i)
      {
         dataOutRef[i] = cooNew[i][0]+cooNew[i][1]*shift+cooNew[i][2]/shift;
      }
   }

   // Transpose op
   auto comm = std::make_shared<Comm<double>>();
   auto transposeOp =
      std::make_unique<Op<View<double, outTy>, View<double, inTy>, p201_t>>(
         comm);

   transposeOp->apply(viewOut, viewIn);

   // check
   if (ranks == 1)
   {
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
   for (std::size_t i = 0; i < dataOutRef.size(); ++i)
   {
      CHECK(dataOut[i] == dataOutRef[i]);
   }
}


TEST_CASE("Mpi S1CLCSC3D to DCCSC3D 210", "MpiS1CLCSC3DtoDCCSC3D210")
{
   int rank, ranks;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &ranks);

   // AL out -> JW in
   // Full data and type
   constexpr size_t M = 4;
   constexpr size_t N = 2;
   constexpr size_t K = 3;

   // view
   constexpr std::uint32_t vRank = 3;
   std::array<std::uint32_t, vRank> dimensionsIn{M, N, K};
   std::array<std::uint32_t, vRank> dimensionsOut{N, K, M};

   ptrAndIdx metaIn;
   ptrAndIdx metaOut;

   if (rank == 0 && ranks == 1)
   {
      // Populate meta for fully populated tensor
      // AL space (Stage::PPM and Stage::MPM)
      metaIn = Index::densePtrAndIdx<DCCSC3D>(dimensionsIn);
      // Populate meta for fully populated tensor
      // Spectral(JW) space (Stage::PMM and Stage::MMM)
      metaOut = Index::densePtrAndIdxStep1<DCCSC3D>(dimensionsOut);

   }

    constexpr size_t S = (M + (M-1) + (M-2)) * N ;
    std::array<double, S> dataIn = {/*k0*/ 1, 2, 3, 4,
        /*k0*/ 5, 6, 7, 8,
        /*k1*/ 9, 10, 11,
        /*k1*/ 12, 13, 14,
        /*k2*/ 15, 16,
        /*k2*/ 17, 18};

    // perm = [2 0 1] -> N K M
    std::array<double, S> dataRef = {/*m0*/ 1, 5,
        /*m0*/ 9,  12,
        /*m0*/ 15, 17,
        /*m1*/ 2,  6,
        /*m1*/ 10, 13,
        /*m1*/ 16, 18,
        /*m2*/ 3,  7,
        /*m2*/ 11, 14,
        /*m3*/ 4,
        /*m3*/ 8};

    std::array<double, S> dataOut;

    // skip setting up pointers and indices
    // they are not used in the serial aka dense tranpose
   std::array<std::vector<std::uint32_t>, vRank> pointersIn = {
      {{}, metaIn.ptr, {}}};
   std::array<std::vector<std::uint32_t>, vRank> indicesIn = {
      {{}, metaIn.idx, {}}};
   std::array<std::vector<std::uint32_t>, vRank> pointersOut = {
      {{}, metaOut.ptr, {}}};
   std::array<std::vector<std::uint32_t>, vRank> indicesOut = {
      {{}, metaOut.idx, {}}};
    using inTy = S1CLCSC3D;
    using outTy = DCCSC3D;
    View<double, inTy> viewIn(dataIn, dimensionsIn, pointersIn, indicesIn);
    View<double, outTy> viewOut(dataOut, dimensionsOut, pointersOut, indicesOut);

    // Transpose op
    using namespace QuICC::Transpose::Mpi;
    using namespace QuICC::Transpose;

    auto comm = std::make_shared<Comm<double>>();
    auto transposeOp =
      std::make_unique<Op<View<double, outTy>, View<double, inTy>, p201_t>>(comm);

    transposeOp->apply(viewOut, viewIn);

    // check
    for (std::uint64_t s = 0; s < S; ++s)
    {
        CHECK(dataRef[s] == viewOut[s]);
    }
}
