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


using namespace QuICC::Transpose::Mpi;


int main(int argc, char** argv)
{
   MPI_Init(NULL, NULL);

   int rank, ranks;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &ranks);

   // if(rank == 0)
   // {
   //     using namespace std::chrono_literals;
   //     volatile bool wait = true;
   //     char hostname[HOST_NAME_MAX];
   //     gethostname(hostname, HOST_NAME_MAX);
   //     std::cerr << "PID " << getpid() << " on " << hostname
   //         << " ready to attach" << std::endl;
   //     while (wait == true)
   //     std::this_thread::sleep_for(1s);
   // }
   // MPI_Barrier(MPI_COMM_WORLD);

   Catch::Session session; // There must be exactly one instance

   auto returnCode = session.run();

   MPI_Finalize();

   return returnCode;
}

using namespace QuICC::Memory;
using namespace QuICC::View;

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
   Index::ptrAndIdx metaIn;
   Index::ptrAndIdx metaOut;

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
