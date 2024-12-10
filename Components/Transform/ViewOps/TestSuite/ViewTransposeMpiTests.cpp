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

#include "Environment/QuICCEnv.hpp"
#include "TestSuite/ViewMeta.hpp"
#include "ViewOps/Transpose/Op.hpp"
#include "ViewOps/ViewIndexUtils.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#include "ViewOps/ViewSizeUtils.hpp"


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

TEST_CASE("Mpi DCCSC3D to DCCSC3D 201", "MpiDCCSC3DtoDCCSC3D201")
{
   int rank, ranks;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &ranks);

   assert(ranks == 1 || ranks == 2 || ranks == 4);

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

   using inTy = DCCSC3D;
   using outTy = DCCSC3D;

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
      metaIn = Index::densePtrAndIdx<inTy>(dimensionsIn);
      // Populate meta for fully populated tensor
      // AL space (Stage::PPM and Stage::MPM)
      metaOut = Index::densePtrAndIdx<outTy>(dimensionsOut);
   }
   else if (rank == 0 && ranks == 2)
   {
      metaIn.ptr = {0, 2, 2};
      metaIn.idx = {0, 1};
      dataIn.resize(metaIn.idx.size() * M);
      // perm = [2 0 1] -> N K M
      metaOut.ptr = {0, 2, 4, 4, 4};
      metaOut.idx = {0, 1, 0, 1};
      dataOut.resize(metaOut.idx.size() * N);
      dataOutRef.resize(metaOut.idx.size() * N);
   }
   else if (rank == 1 && ranks == 2)
   {
      metaIn.ptr = {0, 0, 2};
      metaIn.idx = {0, 1};
      dataIn.resize(metaIn.idx.size() * M);
      // perm = [2 0 1] -> N K M
      metaOut.ptr = {0, 0, 0, 2, 4};
      metaOut.idx = {0, 1, 0, 1};
      dataOut.resize(metaOut.idx.size() * N);
      dataOutRef.resize(metaOut.idx.size() * N);
   }
   else if (ranks == 4)
   {
      std::string path = "_refdata/Framework/LoadSplitter/WLFl/";
      std::string dist = "Tubular";
      std::string id = "103";
      auto setup = readDimsAndMeta(path, dist, id);

      dimensionsIn = {setup.modsDims[2], setup.physDims[1], setup.physDims[0]};
      dimensionsOut = {setup.physDims[1], setup.physDims[0], setup.modsDims[2]};

      metaIn.ptr = setup.metaFT.ptr;
      metaIn.idx = setup.metaFT.idx;
      auto sizeIn = getDataSize<inTy>(dimensionsIn, metaIn);
      dataIn.resize(sizeIn);
      // perm = [2 0 1] MLN -> LNM
      metaOut.ptr = setup.metaAL.ptr;
      metaOut.idx = setup.metaAL.idx;
      auto sizeOut = getDataSize<outTy>(dimensionsOut, metaOut);
      dataOut.resize(sizeOut);
      dataOutRef.resize(sizeOut);
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
         dataIn[i] = cooOld[i][0] + cooOld[i][1] * shift + cooOld[i][2] / shift;
      }
      auto cooNew = ::QuICC::View::getCoo<View<double, outTy>, p201_t>(viewOut);
      for (std::size_t i = 0; i < cooNew.size(); ++i)
      {
         dataOutRef[i] =
            cooNew[i][0] + cooNew[i][1] * shift + cooNew[i][2] / shift;
      }
   }

   // Transpose op
   auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
   auto comm = std::make_shared<Comm<double>>(mem);
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
   else
   {
      for (std::size_t i = 0; i < dataOutRef.size(); ++i)
      {
         CHECK(dataOut[i] == dataOutRef[i]);
      }
   }
}

TEST_CASE("Mpi DCCSC3D to DCCSC3DJIK 201", "MpiDCCSC3DtoDCCSC3DJIK201")
{
   int rank, ranks;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &ranks);

   assert(ranks == 1 || ranks == 2 || ranks == 4);

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

   using inTy = DCCSC3D;
   using outTy = DCCSC3DJIK;

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
      metaIn = Index::densePtrAndIdx<inTy>(dimensionsIn);
      // Populate meta for fully populated tensor
      // AL space (Stage::PPM and Stage::MPM)
      metaOut = Index::densePtrAndIdx<outTy>(dimensionsOut);
   }
   else if (rank == 0 && ranks == 2)
   {
      metaIn.ptr = {0, 2, 2};
      metaIn.idx = {0, 1};
      dataIn.resize(metaIn.idx.size() * M);
      // perm = [2 0 1] -> N K M
      metaOut.ptr = {0, 2, 4, 4, 4};
      metaOut.idx = {0, 1, 0, 1};
      dataOut.resize(metaOut.idx.size() * N);
      dataOutRef.resize(metaOut.idx.size() * N);
   }
   else if (rank == 1 && ranks == 2)
   {
      metaIn.ptr = {0, 0, 2};
      metaIn.idx = {0, 1};
      dataIn.resize(metaIn.idx.size() * M);
      // perm = [2 0 1] -> N K M
      metaOut.ptr = {0, 0, 0, 2, 4};
      metaOut.idx = {0, 1, 0, 1};
      dataOut.resize(metaOut.idx.size() * N);
      dataOutRef.resize(metaOut.idx.size() * N);
   }
   else if (ranks == 4)
   {
      std::string path = "_refdata/Framework/LoadSplitter/WLFl/";
      std::string dist = "Tubular";
      std::string id = "103";
      auto setup = readDimsAndMeta(path, dist, id);

      dimensionsIn = {setup.modsDims[2], setup.physDims[1], setup.physDims[0]};
      dimensionsOut = {setup.physDims[1], setup.physDims[0], setup.modsDims[2]};

      metaIn.ptr = setup.metaFT.ptr;
      metaIn.idx = setup.metaFT.idx;
      auto sizeIn = getDataSize<inTy>(dimensionsIn, metaIn);
      dataIn.resize(sizeIn);
      // perm = [2 0 1] MLN -> LNM
      metaOut.ptr = setup.metaAL.ptr;
      metaOut.idx = setup.metaAL.idx;
      auto sizeOut = getDataSize<outTy>(dimensionsOut, metaOut);
      dataOut.resize(sizeOut);
      dataOutRef.resize(sizeOut);
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
         dataIn[i] = cooOld[i][0] + cooOld[i][1] * shift + cooOld[i][2] / shift;
      }
      auto cooNew = ::QuICC::View::getCoo<View<double, outTy>, p201_t>(viewOut);
      for (std::size_t i = 0; i < cooNew.size(); ++i)
      {
         dataOutRef[i] =
            cooNew[i][0] + cooNew[i][1] * shift + cooNew[i][2] / shift;
      }
   }

   // Transpose op
   auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
   auto comm = std::make_shared<Comm<double>>(mem);
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
               // 201 -> nkm -> JIK -> knm
               auto knm = k + n * K + m * K * N;
               CHECK(viewIn[mnk] == viewOut[knm]);
            }
         }
      }
   }
   else
   {
      for (std::size_t i = 0; i < dataOutRef.size(); ++i)
      {
         CHECK(dataOut[i] == dataOutRef[i]);
      }
   }
}

TEST_CASE("Mpi DCCSC3D to DCCSC3D 120", "MpiDCCSC3DtoDCCSC3D120")
{
   int rank, ranks;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &ranks);

   assert(ranks == 1 || ranks == 4);

   // FFT out -> AL in
   // Full data and type
   constexpr size_t M = 4;
   constexpr size_t N = 2;
   constexpr size_t K = 2;

   // view
   constexpr std::uint32_t vRank = 3;
   std::array<std::uint32_t, vRank> dimensionsIn{M, N, K};
   std::array<std::uint32_t, vRank> dimensionsOut{K, M, N};

   std::vector<double> dataIn;
   std::vector<double> dataOut;
   std::vector<double> dataOutRef;
   ptrAndIdx metaIn;
   ptrAndIdx metaOut;

   using inTy = DCCSC3D;
   using outTy = DCCSC3D;

   // Setup ref data and input data
   using namespace QuICC::Transpose::Mpi;
   using namespace QuICC::Transpose;
   if (rank == 0 && ranks == 1)
   {
      constexpr size_t S = M * N * K;
      dataIn = {/*k0*/ 1, 2, 3, 4,
         /*k0*/ 5, 6, 7, 8,
         /*k1*/ 9, 10, 11, 12,
         /*k1*/ 13, 14, 15, 16};

      // perm = [1 2 0] -> K M N
      dataOut.resize(S);

      // Populate meta for fully populated tensor
      // AL space (Stage::PPM and Stage::MPM)
      metaIn = Index::densePtrAndIdx<inTy>(dimensionsIn);
      // Populate meta for fully populated tensor
      // Physical space (Stage::PPP and Stage::MPP)
      metaOut = Index::densePtrAndIdx<outTy>(dimensionsOut);
   }
   else if (ranks == 4)
   {
      std::string path = "_refdata/Framework/LoadSplitter/WLFl/";
      std::string dist = "Tubular";
      std::string id = "103";
      auto setup = readDimsAndMeta(path, dist, id);

      dimensionsIn = {setup.physDims[1], setup.physDims[0], setup.modsDims[2]};
      dimensionsOut = {setup.modsDims[2], setup.physDims[1], setup.physDims[0]};

      metaIn.ptr = setup.metaAL.ptr;
      metaIn.idx = setup.metaAL.idx;
      auto sizeIn = getDataSize<inTy>(dimensionsIn, metaIn);
      dataIn.resize(sizeIn);
      // perm = [1 2 0] LNM -> MLN
      metaOut.ptr = setup.metaFT.ptr;
      metaOut.idx = setup.metaFT.idx;
      auto sizeOut = getDataSize<outTy>(dimensionsOut, metaOut);
      dataOut.resize(sizeOut);
      dataOutRef.resize(sizeOut);
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
         dataIn[i] = cooOld[i][0] + cooOld[i][1] * shift + cooOld[i][2] / shift;
      }
      auto cooNew = ::QuICC::View::getCoo<View<double, outTy>, p120_t>(viewOut);
      for (std::size_t i = 0; i < cooNew.size(); ++i)
      {
         dataOutRef[i] =
            cooNew[i][0] + cooNew[i][1] * shift + cooNew[i][2] / shift;
      }
   }

   // Transpose op
   auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
   auto comm = std::make_shared<Comm<double>>(mem);
   auto transposeOp =
      std::make_unique<Op<View<double, outTy>, View<double, inTy>, p120_t>>(
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
               auto kmn = k + m * K + n * K * M;
               CHECK(viewIn[mnk] == viewOut[kmn]);
            }
         }
      }
   }
   else
   {
      for (std::size_t i = 0; i < dataOutRef.size(); ++i)
      {
         CHECK(dataOut[i] == dataOutRef[i]);
      }
   }
}

TEST_CASE("Mpi DCCSC3DJIK to DCCSC3D 120", "MpiDCCSC3DJIKtoDCCSC3D120")
{
   int rank, ranks;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &ranks);

   assert(ranks == 1 || ranks == 4);

   // FFT out -> AL in
   // Full data and type
   constexpr size_t M = 4;
   constexpr size_t N = 2;
   constexpr size_t K = 2;

   // view
   constexpr std::uint32_t vRank = 3;
   std::array<std::uint32_t, vRank> dimensionsIn{M, N, K};
   std::array<std::uint32_t, vRank> dimensionsOut{K, M, N};

   std::vector<double> dataIn;
   std::vector<double> dataOut;
   std::vector<double> dataOutRef;
   ptrAndIdx metaIn;
   ptrAndIdx metaOut;

   using inTy = DCCSC3DJIK;
   using outTy = DCCSC3D;

   // Setup ref data and input data
   using namespace QuICC::Transpose::Mpi;
   using namespace QuICC::Transpose;
   if (rank == 0 && ranks == 1)
   {
      constexpr size_t S = M * N * K;
      dataIn = {/*k0*/ 1, 5,
         /*k0*/ 2, 6,
         /*k0*/ 3, 7,
         /*K0*/ 4, 8,
         /*k1*/ 9, 13,
         /*k1*/ 10, 14,
         /*k1*/ 11, 15,
         /*k1*/ 12, 16};

      // perm = [1 2 0] -> K M N
      dataOut.resize(S);

      // Populate meta for fully populated tensor
      // AL space (Stage::PPM and Stage::MPM)
      metaIn = Index::densePtrAndIdx<inTy>(dimensionsIn);
      // Populate meta for fully populated tensor
      // Physical space (Stage::PPP and Stage::MPP)
      metaOut = Index::densePtrAndIdx<outTy>(dimensionsOut);
   }
   else if (ranks == 4)
   {
      std::string path = "_refdata/Framework/LoadSplitter/WLFl/";
      std::string dist = "Tubular";
      std::string id = "103";
      auto setup = readDimsAndMeta(path, dist, id);

      dimensionsIn = {setup.physDims[1], setup.physDims[0], setup.modsDims[2]};
      dimensionsOut = {setup.modsDims[2], setup.physDims[1], setup.physDims[0]};

      metaIn.ptr = setup.metaAL.ptr;
      metaIn.idx = setup.metaAL.idx;
      auto sizeIn = getDataSize<inTy>(dimensionsIn, metaIn);
      dataIn.resize(sizeIn);
      // perm = [1 2 0] LNM -> MLN
      metaOut.ptr = setup.metaFT.ptr;
      metaOut.idx = setup.metaFT.idx;
      auto sizeOut = getDataSize<outTy>(dimensionsOut, metaOut);
      dataOut.resize(sizeOut);
      dataOutRef.resize(sizeOut);
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
         dataIn[i] = cooOld[i][0] + cooOld[i][1] * shift + cooOld[i][2] / shift;
      }
      auto cooNew = ::QuICC::View::getCoo<View<double, outTy>, p120_t>(viewOut);
      for (std::size_t i = 0; i < cooNew.size(); ++i)
      {
         dataOutRef[i] =
            cooNew[i][0] + cooNew[i][1] * shift + cooNew[i][2] / shift;
      }
   }

   // Transpose op
   auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
   auto comm = std::make_shared<Comm<double>>(mem);
   auto transposeOp =
      std::make_unique<Op<View<double, outTy>, View<double, inTy>, p120_t>>(
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
               auto nmk = n + m * N + k * M * N;
               auto kmn = k + m * K + n * K * M;
               CHECK(viewIn[nmk] == viewOut[kmn]);
            }
         }
      }
   }
   else
   {
      for (std::size_t i = 0; i < dataOutRef.size(); ++i)
      {
         CHECK(dataOut[i] == dataOutRef[i]);
      }
   }
}

TEST_CASE("Mpi S1CLCSC3D to DCCSC3D 201", "MpiS1CLCSC3DtoDCCSC3D201")
{
   int rank, ranks;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &ranks);

   assert(ranks == 1 || ranks == 4);

   // AL out -> JW in
   // Full data and type
   constexpr size_t M = 4;
   constexpr size_t N = 2;
   constexpr size_t K = 3;

   // view
   constexpr std::uint32_t vRank = 3;
   std::array<std::uint32_t, vRank> dimensionsIn{M, N, K};
   std::array<std::uint32_t, vRank> dimensionsOut{N, K, M};

   std::vector<double> dataIn;
   std::vector<double> dataOut;
   std::vector<double> dataOutRef;
   ptrAndIdx metaIn;
   ptrAndIdx metaOut;

   using inTy = S1CLCSC3D;
   using outTy = DCCSC3D;

   if (rank == 0 && ranks == 1)
   {
      dataIn = {
         /*k0*/ 1, 2, 3, 4,
         /*k0*/ 5, 6, 7, 8,
         /*k1*/ 9, 10, 11,
         /*k1*/ 12, 13, 14,
         /*k2*/ 15, 16,
         /*k2*/ 17, 18
      };

      // perm = [2 0 1] -> N K M
      dataOutRef = {
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

      // Populate meta for fully populated tensor
      // AL space (Stage::PPM and Stage::MPM)
      metaIn = Index::densePtrAndIdx<inTy>(dimensionsIn);
      // Populate meta for fully populated tensor
      // Spectral(JW) space (Stage::PMM and Stage::MMM)
      metaOut = Index::densePtrAndIdxStep1<outTy>(dimensionsOut);
      dataOut.resize(dataOutRef.size());
   }
   else if (ranks == 4)
   {
      std::string path = "_refdata/Framework/LoadSplitter/WLFl/";
      std::string dist = "Tubular";
      std::string id = "103";
      auto setup = readDimsAndMeta(path, dist, id);

      dimensionsIn = {setup.modsDims[1], setup.physDims[0], setup.modsDims[2]};
      dimensionsOut = {setup.physDims[0], setup.modsDims[2], setup.modsDims[1]};

      metaIn.ptr = setup.metaAL.ptr;
      metaIn.idx = setup.metaAL.idx;
      auto sizeIn = getDataSize<inTy>(dimensionsIn, metaIn);
      dataIn.resize(sizeIn);
      // perm = [2 0 1] MLN -> LNM
      metaOut.ptr = setup.metaJW.ptr;
      metaOut.idx = setup.metaJW.idx;
      auto sizeOut = getDataSize<outTy>(dimensionsOut, metaOut);
      dataOut.resize(sizeOut);
      dataOutRef.resize(sizeOut);
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
         dataIn[i] = cooOld[i][0] + cooOld[i][1] * shift + cooOld[i][2] / shift;
      }
      auto cooNew = ::QuICC::View::getCoo<View<double, outTy>, p201_t>(viewOut);
      for (std::size_t i = 0; i < cooNew.size(); ++i)
      {
         dataOutRef[i] =
            cooNew[i][0] + cooNew[i][1] * shift + cooNew[i][2] / shift;
      }
   }

   // Transpose op
   auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
   auto comm = std::make_shared<Comm<double>>(mem);
   auto transposeOp =
      std::make_unique<Op<View<double, outTy>, View<double, inTy>, p201_t>>(
         comm);

   transposeOp->apply(viewOut, viewIn);

   // check
   for (std::uint64_t s = 0; s < dataOutRef.size(); ++s)
   {
      CHECK(dataOut[s] == dataOutRef[s]);
   }
}

TEST_CASE("Mpi S1CLCSC3DJIK to DCCSC3DJIK 201", "MpiS1CLCSC3DJIKtoDCCSC3D201JIK")
{
   int rank, ranks;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &ranks);

   assert(ranks == 1 || ranks == 4);

   // AL out -> JW in
   // Full data and type
   constexpr size_t M = 4;
   constexpr size_t N = 2;
   constexpr size_t K = 3;

   // view
   constexpr std::uint32_t vRank = 3;
   std::array<std::uint32_t, vRank> dimensionsIn{M, N, K};
   std::array<std::uint32_t, vRank> dimensionsOut{N, K, M};

   std::vector<double> dataIn;
   std::vector<double> dataOut;
   std::vector<double> dataOutRef;
   ptrAndIdx metaIn;
   ptrAndIdx metaOut;

   using inTy = S1CLCSC3DJIK;
   using outTy = DCCSC3DJIK;

   if (rank == 0 && ranks == 1)
   {
      dataIn = {
         /*k0*/ 1, 5,
         /*k0*/ 2, 6,
         /*k0*/ 3, 7,
         /*k0*/ 4, 8,
         /*k1*/ 9, 12,
         /*k1*/ 10, 13,
         /*k1*/ 11, 14,
         /*k2*/ 15, 17,
         /*k2*/ 16, 18
      };

      // perm = [2 0 1] -> N K M
      dataOutRef = {
         /*m0*/ 1,
         /*m0*/ 5,
         /*m1*/ 2, 9,
         /*m1*/ 6, 12,
         /*m2*/ 3, 10, 15,
         /*m2*/ 7, 13, 17,
         /*m3*/ 4, 11, 16,
         /*m3*/ 8, 14, 18
      };

      // Populate meta for fully populated tensor
      // AL space (Stage::PPM and Stage::MPM)
      metaIn = Index::densePtrAndIdx<inTy>(dimensionsIn);
      // Populate meta for fully populated tensor
      // Spectral(JW) space (Stage::PMM and Stage::MMM)
      metaOut = Index::densePtrAndIdxStep1<outTy>(dimensionsOut);
      dataOut.resize(dataOutRef.size());
   }
   else if (ranks == 4)
   {
      std::string path = "_refdata/Framework/LoadSplitter/WLFl/";
      std::string dist = "Tubular";
      std::string id = "103";
      auto setup = readDimsAndMeta(path, dist, id);

      dimensionsIn = {setup.modsDims[1], setup.physDims[0], setup.modsDims[2]};
      dimensionsOut = {setup.physDims[0], setup.modsDims[2], setup.modsDims[1]};

      metaIn.ptr = setup.metaAL.ptr;
      metaIn.idx = setup.metaAL.idx;
      auto sizeIn = getDataSize<inTy>(dimensionsIn, metaIn);
      dataIn.resize(sizeIn);
      // perm = [2 0 1] MLN -> LNM
      metaOut.ptr = setup.metaJW.ptr;
      metaOut.idx = setup.metaJW.idx;
      auto sizeOut = getDataSize<outTy>(dimensionsOut, metaOut);
      dataOut.resize(sizeOut);
      dataOutRef.resize(sizeOut);
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
         dataIn[i] = cooOld[i][0] + cooOld[i][1] * shift + cooOld[i][2] / shift;
      }
      auto cooNew = ::QuICC::View::getCoo<View<double, outTy>, p201_t>(viewOut);
      for (std::size_t i = 0; i < cooNew.size(); ++i)
      {
         dataOutRef[i] =
            cooNew[i][0] + cooNew[i][1] * shift + cooNew[i][2] / shift;
      }
   }

   // Transpose op
   auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
   auto comm = std::make_shared<Comm<double>>(mem);
   auto transposeOp =
      std::make_unique<Op<View<double, outTy>, View<double, inTy>, p201_t>>(
         comm);

   transposeOp->apply(viewOut, viewIn);

   // check
   for (std::uint64_t s = 0; s < dataOutRef.size(); ++s)
   {
      CHECK(dataOut[s] == dataOutRef[s]);
   }
}

TEST_CASE("Mpi DCCSC3D to S1CLCSC3D 120", "MpiDCCSC3DtoS1CLCSC3D120")
{
   int rank, ranks;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &ranks);

   assert(ranks == 1 || ranks == 4);

   // JW out -> AL in
   // Full data and type
   constexpr size_t M = 4;
   constexpr size_t N = 2;
   constexpr size_t K = 3;

   // view
   constexpr std::uint32_t vRank = 3;
   std::array<std::uint32_t, vRank> dimensionsIn{N, K, M};
   std::array<std::uint32_t, vRank> dimensionsOut{M, N, K};

   std::vector<double> dataIn;
   std::vector<double> dataOut;
   std::vector<double> dataOutRef;
   ptrAndIdx metaIn;
   ptrAndIdx metaOut;

   using inTy = DCCSC3D;
   using outTy = S1CLCSC3D;

   if (rank == 0 && ranks == 1)
   {
      // N K M
      dataIn = {
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
      dataOutRef = {
         /*k0*/ 1, 2, 3, 4,
         /*k0*/ 5, 6, 7, 8,
         /*k1*/ 9, 10, 11,
         /*k1*/ 12, 13, 14,
         /*k2*/ 15, 16,
         /*k2*/ 17, 18
      };

      // Populate meta for fully populated tensor
      // Spectral(JW) space (Stage::PMM and Stage::MMM)
      metaIn = Index::densePtrAndIdxStep1<inTy>(dimensionsIn);
      // Populate meta for fully populated tensor
      // AL space (Stage::PPM and Stage::MPM)
      metaOut = Index::densePtrAndIdx<outTy>(dimensionsOut);
      constexpr size_t S = (M + (M - 1) + (M - 2)) * N;

      // dataOut.resize(metaOut.idx.size()*M);
      dataOut.resize(S);
   }
   else if (ranks == 4)
   {
      std::string path = "_refdata/Framework/LoadSplitter/WLFl/";
      std::string dist = "Tubular";
      std::string id = "103";
      auto setup = readDimsAndMeta(path, dist, id);

      dimensionsIn = {setup.physDims[0], setup.modsDims[2], setup.modsDims[1]};
      dimensionsOut = {setup.modsDims[1], setup.physDims[0], setup.modsDims[2]};

      metaIn.ptr = setup.metaJW.ptr;
      metaIn.idx = setup.metaJW.idx;
      auto sizeIn = getDataSize<inTy>(dimensionsIn, metaIn);
      dataIn.resize(sizeIn);
      // perm = [1 2 0] LNM -> MLN
      metaOut.ptr = setup.metaAL.ptr;
      metaOut.idx = setup.metaAL.idx;
      auto sizeOut = getDataSize<outTy>(dimensionsOut, metaOut);
      dataOut.resize(sizeOut);
      dataOutRef.resize(sizeOut);
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

   View<double, inTy> viewIn(dataIn, dimensionsIn, pointersIn, indicesIn);
   View<double, outTy> viewOut(dataOut, dimensionsOut, pointersOut, indicesOut);

   // Transpose op
   using namespace QuICC::Transpose::Mpi;
   using namespace QuICC::Transpose;

   auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
   auto comm = std::make_shared<Comm<double>>(mem);
   auto transposeOp =
      std::make_unique<Op<View<double, outTy>, View<double, inTy>, p120_t>>(
         comm);

   transposeOp->apply(viewOut, viewIn);

   // check
   for (std::uint64_t s = 0; s < dataOutRef.size(); ++s)
   {
      CHECK(dataOut[s] == dataOutRef[s]);
   }
}

TEST_CASE("Mpi DCCSC3DJIK to S1CLCSC3DJIK 120", "MpiDCCSC3DJIKtoS1CLCSC3DJIK120")
{
   int rank, ranks;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &ranks);

   assert(ranks == 1 || ranks == 4);

   // JW out -> AL in
   // Full data and type
   constexpr size_t M = 4;
   constexpr size_t N = 2;
   constexpr size_t K = 3;

   // view
   constexpr std::uint32_t vRank = 3;
   std::array<std::uint32_t, vRank> dimensionsIn{N, K, M};
   std::array<std::uint32_t, vRank> dimensionsOut{M, N, K};

   std::vector<double> dataIn;
   std::vector<double> dataOut;
   std::vector<double> dataOutRef;
   ptrAndIdx metaIn;
   ptrAndIdx metaOut;

   using inTy = DCCSC3DJIK;
   using outTy = S1CLCSC3DJIK;

   if (rank == 0 && ranks == 1)
   {
      // N K M
      dataIn = {
         /*m0*/ 1,
         /*m0*/ 5,
         /*m1*/ 2, 9,
         /*m1*/ 6, 12,
         /*m2*/ 3, 10, 15,
         /*m2*/ 7, 13, 17,
         /*m3*/ 4, 11, 16,
         /*m3*/ 8, 14, 18
      };

      // perm = [1 2 0] -> M N K
      dataOutRef = {
         /*k0*/ 1, 5,
         /*k0*/ 2, 6,
         /*k0*/ 3, 7,
         /*k0*/ 4, 8,
         /*k1*/ 9, 12,
         /*k1*/ 10, 13,
         /*k1*/ 11, 14,
         /*k2*/ 15, 17,
         /*k2*/ 16, 18
      };

      // Populate meta for fully populated tensor
      // Spectral(JW) space (Stage::PMM and Stage::MMM)
      metaIn = Index::densePtrAndIdxStep1<inTy>(dimensionsIn);
      // Populate meta for fully populated tensor
      // AL space (Stage::PPM and Stage::MPM)
      metaOut = Index::densePtrAndIdx<outTy>(dimensionsOut);
      constexpr size_t S = (M + (M - 1) + (M - 2)) * N;

      // dataOut.resize(metaOut.idx.size()*M);
      dataOut.resize(S);
   }
   else if (ranks == 4)
   {
      std::string path = "_refdata/Framework/LoadSplitter/WLFl/";
      std::string dist = "Tubular";
      std::string id = "103";
      auto setup = readDimsAndMeta(path, dist, id);

      dimensionsIn = {setup.physDims[0], setup.modsDims[2], setup.modsDims[1]};
      dimensionsOut = {setup.modsDims[1], setup.physDims[0], setup.modsDims[2]};

      metaIn.ptr = setup.metaJW.ptr;
      metaIn.idx = setup.metaJW.idx;
      auto sizeIn = getDataSize<inTy>(dimensionsIn, metaIn);
      dataIn.resize(sizeIn);
      // perm = [1 2 0] LNM -> MLN
      metaOut.ptr = setup.metaAL.ptr;
      metaOut.idx = setup.metaAL.idx;
      auto sizeOut = getDataSize<outTy>(dimensionsOut, metaOut);
      dataOut.resize(sizeOut);
      dataOutRef.resize(sizeOut);
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

   View<double, inTy> viewIn(dataIn, dimensionsIn, pointersIn, indicesIn);
   View<double, outTy> viewOut(dataOut, dimensionsOut, pointersOut, indicesOut);

   // Transpose op
   using namespace QuICC::Transpose::Mpi;
   using namespace QuICC::Transpose;

   auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
   auto comm = std::make_shared<Comm<double>>(mem);
   auto transposeOp =
      std::make_unique<Op<View<double, outTy>, View<double, inTy>, p120_t>>(
         comm);

   transposeOp->apply(viewOut, viewIn);

   // check
   for (std::uint64_t s = 0; s < dataOutRef.size(); ++s)
   {
      CHECK(dataOut[s] == dataOutRef[s]);
   }
}

