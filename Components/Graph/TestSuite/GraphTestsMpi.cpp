#define CATCH_CONFIG_RUNNER

#include <catch2/catch.hpp>
#include <memory>

// QuICC
#include "Environment/QuICCEnv.hpp"
#include "Graph/Shims/MlirShims.hpp"
#include "Graph/BackendsMap.hpp"
#include "Graph/OpsMap.hpp"
#include "Graph/Jit.hpp"
#include "Memory/Cpu/NewDelete.hpp"
#include "Memory/Cuda/Malloc.hpp"
#include "Memory/Memory.hpp"

#include "Types/Math.hpp"
#include "Profiler/Interface.hpp"
#include "TestSuite/Io.hpp"
#include "Utils.hpp"

int main(int argc, char **argv)
{
  QuICC::QuICCEnv();

  QuICC::Profiler::Initialize();

  Catch::Session session; // There must be exactly one instance

  auto returnCode = session.run();

  QuICC::Profiler::Finalize();

  return returnCode;
}


TEST_CASE("Parallel 3D Fwd", "[Parallel3DFwd]")
{
  std::string inputFilename = "./simple-3d-fwd.mlir";
  llvm::ErrorOr<std::unique_ptr<llvm::MemoryBuffer>> fileOrErr =
      llvm::MemoryBuffer::getFileOrSTDIN(inputFilename);
  if (std::error_code EC = fileOrErr.getError()) {
    llvm::errs() << "Could not open input file: " << EC.message() << "\n";
    CHECK(false);
  }

  // Parse the input mlir.
  llvm::SourceMgr sourceMgr;
  sourceMgr.AddNewSourceBuffer(std::move(*fileOrErr), llvm::SMLoc());

  // Load size meta info
  #if 0
  std::string dist = "Serial";
  int rank = 0;
  int ranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  assert(rank == 0);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  assert(ranks == 1);
  #else
  std::string dist = "Tubular";
  int rank = 0;
  int ranks = 4;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  assert(ranks == 4);
  #endif
  std::string path = "/home/gcastigl/codes/QuICC/build/Components/Framework/TestSuite/_refdata/Framework/LoadSplitter/WLFl/";
  std::string fileName = path+dist+"/P_id103_np"+std::to_string(ranks)+"_r"+std::to_string(rank)+"_stage0_meta.dat";
  std::vector<double> dbJW;
  QuICC::TestSuite::readList(dbJW, fileName);
  std::uint32_t modsK = dbJW[0];
  std::uint32_t K = dbJW[1];

  fileName = path+dist+"/P_id103_np"+std::to_string(ranks)+"_r"+std::to_string(rank)+"_stage1_meta.dat";
   std::vector<double> dbAL;
  QuICC::TestSuite::readList(dbAL, fileName);
  std::uint32_t modsN = dbAL[0];
  std::uint32_t N = dbAL[1];

  fileName = path+dist+"/P_id103_np"+std::to_string(ranks)+"_r"+std::to_string(rank)+"_stage2_meta.dat";
   std::vector<double> dbFr;
  QuICC::TestSuite::readList(dbFr, fileName);
  std::uint32_t modsM = dbFr[0];
  std::uint32_t M = dbFr[1];

  // Grid dimensions
  constexpr std::uint32_t dim = 3;
  std::array<std::uint32_t, dim> physDims{M, N, K};
  std::array<std::uint32_t, dim> modsDims{modsM, modsN, modsK};

  // Unpack metadata Fourier stage
  using namespace QuICC::Graph;
  auto metaFr = unPackMeta(dbFr, K);
  if (rank==0) {
    for (std::size_t i = 0; i < metaFr.ptr.size();++i) {
      std::cout << metaFr.ptr[i] << '\n';
    }
  }

  // Unpack metadata AL stage
  auto metaAL = unPackMeta(dbAL, modsM);

  // Unpack metadata JW stage
  auto metaJW = unPackMetaJW(dbJW, modsN);

  std::array<std::array<std::string, 2>, 3> layOpt;
  layOpt[0] = {"DCCSC3D", "DCCSC3D"};
  layOpt[1] = {"DCCSC3D", "S1CLCSC3D"};
  layOpt[2] = {"DCCSC3D", "DCCSC3D"};

  // Physical space (Stage::PPP and Stage::MPP)
  std::array<std::vector<std::uint32_t>, dim> pointersPhys {{{}, metaFr.ptr, {}}};
  std::array<std::vector<std::uint32_t>, dim> indicesPhys {{{}, metaFr.idx, {}}};

  // Modal space (aka JW space, Stage::PMM and Stage::MMM, QuICC Stage0)
  std::array<std::vector<std::uint32_t>, dim> pointersMods {{{}, metaJW.ptr, {}}};
  std::array<std::vector<std::uint32_t>, dim> indicesMods {{{}, metaJW.idx, {}}};

  // Store meta stages to pass to Jitter
  std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
  meta.push_back({metaFr.ptr.data(), metaFr.ptr.size()});
  meta.push_back({metaFr.idx.data(), metaFr.idx.size()});
  meta.push_back({metaAL.ptr.data(), metaAL.ptr.size()});
  meta.push_back({metaAL.idx.data(), metaAL.idx.size()});
  meta.push_back({metaJW.ptr.data(), metaJW.ptr.size()});
  meta.push_back({metaJW.idx.data(), metaJW.idx.size()});


  QuICC::QuICCEnv().synchronize();
  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  using namespace QuICC::Graph;
  Jit<dim> Jitter(std::move(sourceMgr), mem, physDims, modsDims, layOpt, Stage::MMM, Stage::PPP, meta);

  // host mem block
  std::size_t physS = M*indicesPhys[1].size();
  QuICC::Memory::MemBlock<double> R(physS, mem.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOut(modsK*metaJW.idx.size(), mem.get());

  // host view
  R_DCCSC3D_t RView({R.data(), R.size()}, physDims, pointersPhys, indicesPhys);
  C_DCCSC3D_t modsOutView({modsOut.data(), modsOut.size()}, {modsK, modsM, modsN}, pointersMods, indicesMods);

  // set input phys
  double val = 1.0;
  for(std::size_t m = 0; m < RView.size(); ++m)
  {
    RView[m] = val;
  }

  // Apply graph
  Jitter.apply(modsOutView, RView);

  // Check
  double eps = 1e-14;
  CHECK(std::abs(modsOutView[0].real() - sqrt(2.0)*QuICC::Math::PI) <= eps);
  CHECK(std::abs(modsOutView[0].imag()) <= eps);
  for(std::size_t m = 1; m < modsOutView.size(); ++m)
  {
    CHECK(std::abs(modsOutView[m]) <= eps);
  }

}

// TEST_CASE("Serial 3D Bwd", "[Serial3DBwd]")
// {
//   std::string inputFilename = "./simple-3d-bwd.mlir";
//   llvm::ErrorOr<std::unique_ptr<llvm::MemoryBuffer>> fileOrErr =
//       llvm::MemoryBuffer::getFileOrSTDIN(inputFilename);
//   if (std::error_code EC = fileOrErr.getError()) {
//     llvm::errs() << "Could not open input file: " << EC.message() << "\n";
//     CHECK(false);
//   }

//   // Parse the input mlir.
//   llvm::SourceMgr sourceMgr;
//   sourceMgr.AddNewSourceBuffer(std::move(*fileOrErr), llvm::SMLoc());

//   // Load size meta info
//   #if 1
//   std::string dist = "Serial";
//   int ranks = 1;
//   #else
//   std::string dist = "Tubular";
//   int ranks = 4;
//   #endif
//   std::string path = "/home/gcastigl/codes/QuICC/build/Components/Framework/TestSuite/_data/Framework/LoadSplitter/WLFl/";
//   std::string fileName = path+dist+"/P_id103_np"+std::to_string(ranks)+"_r0_stage0_meta.dat";
//   std::vector<double> dbJW;
//   QuICC::TestSuite::readList(dbJW, fileName);
//   std::uint32_t modsK = dbJW[0];
//   std::uint32_t K = dbJW[1];

//   fileName = path+dist+"/P_id103_np"+std::to_string(ranks)+"_r0_stage1_meta.dat";
//    std::vector<double> dbAL;
//   QuICC::TestSuite::readList(dbAL, fileName);
//   std::uint32_t modsN = dbAL[0];
//   std::uint32_t N = dbAL[1];

//   fileName = path+dist+"/P_id103_np"+std::to_string(ranks)+"_r0_stage2_meta.dat";
//    std::vector<double> dbFr;
//   QuICC::TestSuite::readList(dbFr, fileName);
//   std::uint32_t modsM = dbFr[0];
//   std::uint32_t M = dbFr[1];

//   // Grid dimensions
//   constexpr std::uint32_t rank = 3;
//   std::array<std::uint32_t, rank> physDims{M, N, K};
//   std::array<std::uint32_t, rank> modsDims{modsM, modsN, modsK};

//   // Unpack metadata Fourier stage
//   using namespace QuICC::Graph;
//   auto metaFr = unPackMeta(dbFr, K);

//   // Unpack metadata AL stage
//   auto metaAL = unPackMeta(dbAL, modsM);

//   // Unpack metadata JW stage
//   auto metaJW = unPackMetaJW(dbJW, modsN);

//   std::array<std::array<std::string, 2>, 3> layOpt;
//   layOpt[0] = {"DCCSC3D", "DCCSC3D"};
//   layOpt[1] = {"DCCSC3D", "S1CLCSC3D"};
//   layOpt[2] = {"DCCSC3D", "DCCSC3D"};

//   // Physical space
//   std::array<std::vector<std::uint32_t>, rank> pointersPhys {{{}, metaFr.ptr, {}}};
//   std::array<std::vector<std::uint32_t>, rank> indicesPhys {{{}, metaFr.idx, {}}};

//   // Modal space
//   std::array<std::vector<std::uint32_t>, rank> pointersMods {{{}, metaJW.ptr, {}}};
//   std::array<std::vector<std::uint32_t>, rank> indicesMods {{{}, metaJW.idx, {}}};

//   // Store meta stages to pass to Jitter
//   std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
//   meta.push_back({metaFr.ptr.data(), metaFr.ptr.size()});
//   meta.push_back({metaFr.idx.data(), metaFr.idx.size()});
//   meta.push_back({metaAL.ptr.data(), metaAL.ptr.size()});
//   meta.push_back({metaAL.idx.data(), metaAL.idx.size()});
//   meta.push_back({metaJW.ptr.data(), metaJW.ptr.size()});
//   meta.push_back({metaJW.idx.data(), metaJW.idx.size()});

//   auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
//   using namespace QuICC::Graph;
//   Jit<rank> Jitter(std::move(sourceMgr), mem, physDims, modsDims, layOpt, Stage::PPP, Stage::MMM, meta);

//   // host mem block
//   std::size_t physS = M*indicesPhys[1].size();
//   QuICC::Memory::MemBlock<double> R(physS, mem.get());
//   QuICC::Memory::MemBlock<std::complex<double>> modsIn(modsK*indicesMods[1].size(), mem.get());

//   // host view
//   R_DCCSC3D_t RView({R.data(), R.size()}, physDims, pointersPhys, indicesPhys);
//   C_DCCSC3D_t modsInView({modsIn.data(), modsIn.size()}, {modsK, modsM, modsN}, pointersMods, indicesMods);

//   // set input modes
//   std::complex<double> val = {sqrt(2.0)/2.0, 0.0};
//   modsInView[0] = val;
//   for(std::size_t m = 1; m < modsInView.size(); ++m)
//   {
//     modsInView[m] = {0.0, 0.0};
//   }

//   // Apply graph
//   Jitter.apply(RView, modsInView);

//   // Check
//   double eps = 1e-15;
//   for(std::size_t m = 0; m < RView.size(); ++m)
//   {
//     CHECK(std::abs(RView[m] - 0.5/QuICC::Math::PI) <= eps);
//   }

// }

// TEST_CASE("Serial 3D Loop", "[Serial3DLoop]")
// {
//   std::string inputFilename = "./simple-3d-loop.mlir";
//   llvm::ErrorOr<std::unique_ptr<llvm::MemoryBuffer>> fileOrErr =
//       llvm::MemoryBuffer::getFileOrSTDIN(inputFilename);
//   if (std::error_code EC = fileOrErr.getError()) {
//     llvm::errs() << "Could not open input file: " << EC.message() << "\n";
//     CHECK(false);
//   }

//   // Parse the input mlir.
//   llvm::SourceMgr sourceMgr;
//   sourceMgr.AddNewSourceBuffer(std::move(*fileOrErr), llvm::SMLoc());

//   // Load size meta info
//   #if 1
//   std::string dist = "Serial";
//   int ranks = 1;
//   #else
//   std::string dist = "Tubular";
//   int ranks = 4;
//   #endif
//   std::string path = "/home/gcastigl/codes/QuICC/build/Components/Framework/TestSuite/_data/Framework/LoadSplitter/WLFl/";
//   std::string fileName = path+dist+"/P_id103_np"+std::to_string(ranks)+"_r0_stage0_meta.dat";
//   std::vector<double> dbJW;
//   QuICC::TestSuite::readList(dbJW, fileName);
//   std::uint32_t modsK = dbJW[0];
//   std::uint32_t K = dbJW[1];

//   fileName = path+dist+"/P_id103_np"+std::to_string(ranks)+"_r0_stage1_meta.dat";
//    std::vector<double> dbAL;
//   QuICC::TestSuite::readList(dbAL, fileName);
//   std::uint32_t modsN = dbAL[0];
//   std::uint32_t N = dbAL[1];

//   fileName = path+dist+"/P_id103_np"+std::to_string(ranks)+"_r0_stage2_meta.dat";
//    std::vector<double> dbFr;
//   QuICC::TestSuite::readList(dbFr, fileName);
//   std::uint32_t modsM = dbFr[0];
//   std::uint32_t M = dbFr[1];

//   // Grid dimensions
//   constexpr std::uint32_t rank = 3;
//   std::array<std::uint32_t, rank> physDims{M, N, K};
//   std::array<std::uint32_t, rank> modsDims{modsM, modsN, modsK};

//   // Unpack metadata Fourier stage
//   using namespace QuICC::Graph;
//   auto metaFr = unPackMeta(dbFr, K);

//   // Unpack metadata AL stage
//   auto metaAL = unPackMeta(dbAL, modsM);

//   // Unpack metadata JW stage
//   auto metaJW = unPackMetaJW(dbJW, modsN);

//   std::array<std::array<std::string, 2>, 3> layOpt;
//   layOpt[0] = {"DCCSC3D", "DCCSC3D"};
//   layOpt[1] = {"DCCSC3D", "S1CLCSC3D"};
//   layOpt[2] = {"DCCSC3D", "DCCSC3D"};

//   // Physical space
//   std::array<std::vector<std::uint32_t>, rank> pointersPhys {{{}, metaFr.ptr, {}}};
//   std::array<std::vector<std::uint32_t>, rank> indicesPhys {{{}, metaFr.idx, {}}};

//   // Modal space
//   std::array<std::vector<std::uint32_t>, rank> pointersMods {{{}, metaJW.ptr, {}}};
//   std::array<std::vector<std::uint32_t>, rank> indicesMods {{{}, metaJW.idx, {}}};

//   // Store meta stages to pass to Jitter
//   std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
//   meta.push_back({metaFr.ptr.data(), metaFr.ptr.size()});
//   meta.push_back({metaFr.idx.data(), metaFr.idx.size()});
//   meta.push_back({metaAL.ptr.data(), metaAL.ptr.size()});
//   meta.push_back({metaAL.idx.data(), metaAL.idx.size()});
//   meta.push_back({metaJW.ptr.data(), metaJW.ptr.size()});
//   meta.push_back({metaJW.idx.data(), metaJW.idx.size()});

//   auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
//   using namespace QuICC::Graph;
//   Jit<rank> Jitter(std::move(sourceMgr), mem, physDims, modsDims, layOpt, Stage::MMM, Stage::MMM, meta);

//   // host mem block
//   QuICC::Memory::MemBlock<std::complex<double>> modsIn(modsK*indicesMods[1].size(), mem.get());
//   QuICC::Memory::MemBlock<std::complex<double>> modsOut(modsK*indicesMods[1].size(), mem.get());

//   // host view
//   C_DCCSC3D_t modsInView({modsIn.data(), modsIn.size()}, {modsK, modsM, modsN}, pointersMods, indicesMods);
//   C_DCCSC3D_t modsOutView({modsOut.data(), modsOut.size()}, {modsK, modsM, modsN}, pointersMods, indicesMods);

//   // set input modes
//   std::complex<double> val = {1.0, 0.0};

//   for(std::size_t m = 0; m < modsInView.size(); ++m)
//   {
//     modsInView[m] = {0.0, 0.0};
//   }
//   modsInView(0,0,0) = val;
//   modsInView(1,0,0) = val;

//   // Apply graph
//   Jitter.apply(modsOutView, modsInView);

//   // Check
//   double eps = 1e-15;
//   for(std::size_t m = 0; m < modsOutView.size(); ++m)
//   {
//     CHECK(std::abs(modsOutView[m] - modsInView[m]) <= eps);
//   }

// }

// TEST_CASE("Serial Multi Var 3D Fwd", "[SerialMultiVar3DFwd]")
// {
//   std::string inputFilename = "./complex-3d-fwd.mlir";
//   llvm::ErrorOr<std::unique_ptr<llvm::MemoryBuffer>> fileOrErr =
//       llvm::MemoryBuffer::getFileOrSTDIN(inputFilename);
//   if (std::error_code EC = fileOrErr.getError()) {
//     llvm::errs() << "Could not open input file: " << EC.message() << "\n";
//     CHECK(false);
//   }

//   // Parse the input mlir.
//   llvm::SourceMgr sourceMgr;
//   sourceMgr.AddNewSourceBuffer(std::move(*fileOrErr), llvm::SMLoc());

//   // Grid dimensions
//   constexpr std::uint32_t rank = 3;
//   std::array<std::uint32_t, rank> physDims{10, 6, 3};
//   std::array<std::uint32_t, rank> modsDims{6, 6, 2};

//   std::array<std::array<std::string, 2>, 3> layOpt;
//   layOpt[0] = {"R_DCCSC3D_t", "C_DCCSC3D_t"};
//   layOpt[1] = {"C_DCCSC3D_t", "C_S1CLCSC3D_t"};
//   layOpt[2] = {"C_DCCSC3D_t", "C_DCCSC3D_t"};

//   auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
//   using namespace QuICC::Graph;
//   Jit<rank> Jitter(std::move(sourceMgr), mem, physDims, modsDims, layOpt, Stage::MMM, Stage::PPP);

//   // setup metadata
//   auto M = physDims[0];
//   auto N = physDims[1];
//   auto K = physDims[2];
//   auto modsM = modsDims[0];
//   auto modsN = modsDims[1];
//   auto modsK = modsDims[2];

//   // Populate meta for fully populated tensor
//   // Physical space
//   std::vector<std::uint32_t> ptrPhys(K+1);
//   std::vector<std::uint32_t> idxPhys(K*N);
//   ptrPhys[0] = 0;
//   for (std::size_t i = 1; i < ptrPhys.size(); ++i) {
//       ptrPhys[i] = ptrPhys[i-1]+N;
//   }
//   for (std::size_t i = 0; i < idxPhys.size(); ++i) {
//       idxPhys[i] = i % N;
//   }
//   std::array<std::vector<std::uint32_t>, rank> pointersPhys {{{}, ptrPhys, {}}};
//   std::array<std::vector<std::uint32_t>, rank> indicesPhys {{{}, idxPhys, {}}};

//   // Populate meta for fully populated tensor
//   // Modal space
//   using namespace QuICC::Graph;
//   auto metaMods = denseTransposePtrAndIdx<C_DCCSC3D_t, C_S1CLCSC3D_t>({modsK, modsM, modsN});
//   std::array<std::vector<std::uint32_t>, rank> pointersMods {{{}, metaMods.ptr, {}}};
//   std::array<std::vector<std::uint32_t>, rank> indicesMods {{{}, metaMods.idx, {}}};

//   // host mem block
//   std::size_t physS = M*indicesPhys[1].size();
//   QuICC::Memory::MemBlock<double> R(physS, mem.get());
//   QuICC::Memory::MemBlock<double> Phi(physS, mem.get());
//   QuICC::Memory::MemBlock<double> Theta(physS, mem.get());
//   QuICC::Memory::MemBlock<std::complex<double>> modsOut(modsK*metaMods.idx.size(), mem.get());

//   // host view
//   R_DCCSC3D_t RView({R.data(), R.size()}, physDims, pointersPhys, indicesPhys);
//   R_DCCSC3D_t PhiView({Phi.data(), Phi.size()}, physDims, pointersPhys, indicesPhys);
//   R_DCCSC3D_t ThetaView({Theta.data(), Theta.size()}, physDims, pointersPhys, indicesPhys);
//   C_DCCSC3D_t modsOutView({modsOut.data(), modsOut.size()}, {modsK, modsM, modsN}, pointersMods, indicesMods);

//   // set input phys
//   double val = 1.0;
//   for(std::size_t m = 0; m < RView.size(); ++m)
//   {
//     RView[m] = val;
//     PhiView[m] = val;
//     ThetaView[m] = val;
//   }

//   // Apply graph
//   Jitter.apply(modsOutView, RView, PhiView, ThetaView);

//   // Check
//   double eps = 1e-15;
//   CHECK(std::abs(modsOutView[0].real() - sqrt(2.0)*QuICC::Math::PI) <= eps);
//   CHECK(std::abs(modsOutView[0].imag()) <= eps);
//   for(std::size_t m = 1; m < modsOutView.size(); ++m)
//   {
//     CHECK(std::abs(modsOutView[m]) <= eps);
//   }
// }
