
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


using namespace QuICC::Memory;
using namespace QuICC::View;
using namespace QuICC::TestSuite;


TEST_CASE("Mpi DCCSC3DJIK to S1CLCSC3DJIK 120 Cuda", "MpiDCCSC3DJIKtoS1CLCSC3DJIK120Cuda")
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

   // device mem
   auto memDev = std::make_shared<QuICC::Memory::Cuda::Malloc>();
   QuICC::Memory::MemBlock<double> memBlockIn(dataIn.size(), memDev.get());
   QuICC::Memory::MemBlock<double> memBlockOut(dataOut.size(), memDev.get());

   QuICC::Memory::MemBlock<std::uint32_t> memBlockPtrIn(pointersIn[1].size(), memDev.get());
   QuICC::Memory::MemBlock<std::uint32_t> memBlockIdxIn(indicesIn[1].size(), memDev.get());
   QuICC::Memory::MemBlock<std::uint32_t> memBlockPtrOut(pointersOut[1].size(), memDev.get());
   QuICC::Memory::MemBlock<std::uint32_t> memBlockIdxOut(indicesOut[1].size(), memDev.get());

   // set device pointers and indices
   constexpr std::uint32_t dims = 3;
   ViewBase<std::uint32_t> pointersInDev[dims];
   ViewBase<std::uint32_t> indicesInDev[dims];
   ViewBase<std::uint32_t> pointersOutDev[dims];
   ViewBase<std::uint32_t> indicesOutDev[dims];

   pointersInDev[1] = ViewBase<std::uint32_t>(memBlockPtrIn.data(), memBlockPtrIn.size());
   indicesInDev[1] = ViewBase<std::uint32_t>(memBlockIdxIn.data(), memBlockIdxIn.size());
   pointersOutDev[1] = ViewBase<std::uint32_t>(memBlockPtrOut.data(), memBlockPtrOut.size());
   indicesOutDev[1] = ViewBase<std::uint32_t>(memBlockIdxOut.data(), memBlockIdxOut.size());

   // set device views
   View<double, inTy> viewInDev(memBlockIn.data(), memBlockIn.size(),
      dimensionsIn.data(), pointersInDev, indicesInDev);
   View<double, outTy> viewOutDev(memBlockOut.data(), memBlockOut.size(),
      dimensionsOut.data(), pointersOutDev, indicesOutDev);

   // cpu -> gpu
   cudaErrChk(cudaMemcpy(viewInDev.data(), viewIn.data(), viewIn.size() * sizeof(double),
      cudaMemcpyHostToDevice));
   cudaErrChk(cudaMemcpy(pointersInDev[1].data(), pointersIn[1].data(), pointersIn[1].size() * sizeof(std::uint32_t),
      cudaMemcpyHostToDevice));
   cudaErrChk(cudaMemcpy(indicesInDev[1].data(), indicesIn[1].data(), indicesIn[1].size() * sizeof(std::uint32_t),
      cudaMemcpyHostToDevice));
   cudaErrChk(cudaMemcpy(pointersOutDev[1].data(), pointersOut[1].data(), pointersOut[1].size() * sizeof(std::uint32_t),
      cudaMemcpyHostToDevice));
   cudaErrChk(cudaMemcpy(indicesOutDev[1].data(), indicesOut[1].data(), indicesOut[1].size() * sizeof(std::uint32_t),
      cudaMemcpyHostToDevice));

   // Transpose op
   using namespace QuICC::Transpose::Mpi;
   using namespace QuICC::Transpose;

   auto comm = std::make_shared<Comm<double>>();
   auto transposeOp =
      std::make_unique<Op<View<double, outTy>, View<double, inTy>, p120_t>>(
         comm);

   transposeOp->apply(viewOutDev, viewInDev);

   // gpu -> cpu
   cudaErrChk(cudaMemcpy(viewOut.data(), viewOutDev.data(),
      viewOut.size() * sizeof(double), cudaMemcpyDeviceToHost));

   // check
   for (std::uint64_t s = 0; s < dataOutRef.size(); ++s)
   {
      CHECK(dataOut[s] == dataOutRef[s]);
   }
}

