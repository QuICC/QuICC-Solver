/**
 * @file GraphTestsMpi.hpp
 * @brief mpi/distributed graph test with cpu operators
 */
#define CATCH_CONFIG_RUNNER

#include <catch2/catch.hpp>
#include <memory>

// QuICC
#include "Environment/QuICCEnv.hpp"
#include "Graph/Shims/MlirShims.hpp"
#include "Graph/Jit.hpp"
#include "Memory/Cpu/NewDelete.hpp"
#include "Memory/Cuda/Malloc.hpp"
#include "Memory/Memory.hpp"

#include "Types/Math.hpp"
#include "Profiler/Interface.hpp"
#include "TestSuite/Io.hpp"
#include "TestSuite/ViewMeta.hpp"

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
  int rank = 0;
  int ranks = 1;
  #ifdef QUICC_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  #endif

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
  std::string dist;
  if (rank == 0 && ranks == 1)
  {
    dist = "Serial";
  }
  else if (ranks == 4)
  {
    dist = "Tubular";
  }
  else
  {
    assert(false && "missing setup");
  }
  std::string path = "_refdata/Framework/LoadSplitter/WLFl/";
  std::string id = "103";
  using namespace QuICC::TestSuite;
  auto setup = readDimsAndMeta(path, dist, id);

  // Layouts
  std::array<std::array<std::string, 2>, 3> layOpt;
  layOpt[0] = {"DCCSC3D", "DCCSC3D"};
  layOpt[1] = {"DCCSC3D", "S1CLCSC3D"};
  layOpt[2] = {"DCCSC3D", "DCCSC3D"};

  // Physical space (Stage::PPP and Stage::MPP)
  constexpr std::uint32_t dim = 3;
  std::array<std::vector<std::uint32_t>, dim> pointersPhys {{{}, setup.metaFT.ptr, {}}};
  std::array<std::vector<std::uint32_t>, dim> indicesPhys {{{}, setup.metaFT.idx, {}}};

  // Modal space (aka JW space, Stage::PMM and Stage::MMM, QuICC Stage0)
  std::array<std::vector<std::uint32_t>, dim> pointersMods {{{}, setup.metaJW.ptr, {}}};
  std::array<std::vector<std::uint32_t>, dim> indicesMods {{{}, setup.metaJW.idx, {}}};

  // Store meta stages to pass to Jitter
  std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
  meta.push_back({setup.metaFT.ptr.data(), setup.metaFT.ptr.size()});
  meta.push_back({setup.metaFT.idx.data(), setup.metaFT.idx.size()});
  meta.push_back({setup.metaAL.ptr.data(), setup.metaAL.ptr.size()});
  meta.push_back({setup.metaAL.idx.data(), setup.metaAL.idx.size()});
  meta.push_back({setup.metaJW.ptr.data(), setup.metaJW.ptr.size()});
  meta.push_back({setup.metaJW.idx.data(), setup.metaJW.idx.size()});

  // QuICC::QuICCEnv().synchronize();
  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  using namespace QuICC::Graph;
  Jit<dim> Jitter(std::move(sourceMgr), mem, setup.physDims, setup.modsDims, layOpt, Stage::MMM, Stage::PPP, meta);

  // host mem block
  std::size_t physS = setup.physDims[2]*indicesPhys[1].size();
  QuICC::Memory::MemBlock<double> R(physS, mem.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOut(setup.modsDims[0]*indicesMods[1].size(), mem.get());

  // host view
  std::array<std::uint32_t, dim> physDims{setup.physDims[2], setup.physDims[1], setup.physDims[0]};
  R_DCCSC3D_t RView({R.data(), R.size()}, physDims, pointersPhys, indicesPhys);
  C_DCCSC3D_t modsOutView({modsOut.data(), modsOut.size()}, {setup.modsDims[0], setup.modsDims[2], setup.modsDims[1]}, pointersMods, indicesMods);

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
  std::size_t m = 0;
  try
  {
    // we don't know which rank has mode zero
    modsOutView(0,0,0);
    // if we didn't throw this rank has mode zero
    CHECK(std::abs(modsOutView[m].real() - sqrt(2.0)*QuICC::Math::PI) <= eps);
    CHECK(std::abs(modsOutView[m].imag()) <= eps);
    ++m;
  }
  catch(std::logic_error& e)
  {
    // do nothing
  }

  for(; m < modsOutView.size(); ++m)
  {
    CHECK(std::abs(modsOutView[m]) <= eps);
  }

}

TEST_CASE("Parallel 3D Bwd", "[Parallel3DBwd]")
{
  int rank = 0;
  int ranks = 1;
  #ifdef QUICC_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  #endif

  std::string inputFilename = "./simple-3d-bwd.mlir";
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
  std::string dist;
  if (rank == 0 && ranks == 1)
  {
    dist = "Serial";
  }
  else if (ranks == 4)
  {
    dist = "Tubular";
  }
  else
  {
    assert(false && "missing setup");
  }
  std::string path = "_refdata/Framework/LoadSplitter/WLFl/";
  std::string id = "103";
  using namespace QuICC::TestSuite;
  auto setup = readDimsAndMeta(path, dist, id);

  // Layouts
  std::array<std::array<std::string, 2>, 3> layOpt;
  layOpt[0] = {"DCCSC3D", "DCCSC3D"};
  layOpt[1] = {"DCCSC3D", "S1CLCSC3D"};
  layOpt[2] = {"DCCSC3D", "DCCSC3D"};

  // Physical space (Stage::PPP and Stage::MPP)
  constexpr std::uint32_t dim = 3;
  std::array<std::vector<std::uint32_t>, dim> pointersPhys {{{}, setup.metaFT.ptr, {}}};
  std::array<std::vector<std::uint32_t>, dim> indicesPhys {{{}, setup.metaFT.idx, {}}};

  // Modal space (aka JW space, Stage::PMM and Stage::MMM, QuICC Stage0)
  std::array<std::vector<std::uint32_t>, dim> pointersMods {{{}, setup.metaJW.ptr, {}}};
  std::array<std::vector<std::uint32_t>, dim> indicesMods {{{}, setup.metaJW.idx, {}}};

  // Store meta stages to pass to Jitter
  std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
  meta.push_back({setup.metaFT.ptr.data(), setup.metaFT.ptr.size()});
  meta.push_back({setup.metaFT.idx.data(), setup.metaFT.idx.size()});
  meta.push_back({setup.metaAL.ptr.data(), setup.metaAL.ptr.size()});
  meta.push_back({setup.metaAL.idx.data(), setup.metaAL.idx.size()});
  meta.push_back({setup.metaJW.ptr.data(), setup.metaJW.ptr.size()});
  meta.push_back({setup.metaJW.idx.data(), setup.metaJW.idx.size()});

  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  using namespace QuICC::Graph;
  Jit<dim> Jitter(std::move(sourceMgr), mem, setup.physDims, setup.modsDims, layOpt, Stage::PPP, Stage::MMM, meta);

  // host mem block
  std::size_t physS = setup.physDims[2]*indicesPhys[1].size();
  QuICC::Memory::MemBlock<double> R(physS, mem.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsIn(setup.modsDims[0]*indicesMods[1].size(), mem.get());

  // host view
  std::array<std::uint32_t, dim> physDims{setup.physDims[2], setup.physDims[1], setup.physDims[0]};
  R_DCCSC3D_t RView({R.data(), R.size()}, physDims, pointersPhys, indicesPhys);
  C_DCCSC3D_t modsInView({modsIn.data(), modsIn.size()}, {setup.modsDims[0], setup.modsDims[2], setup.modsDims[1]}, pointersMods, indicesMods);

  // set input modes
  std::complex<double> val = {sqrt(2.0)/2.0, 0.0};
  std::size_t m = 0;
  try
  {
    // we don't know which rank has mode zero
    modsInView(0,0,0);
    // if we didn't throw this rank has mode zero
    modsInView[m] = val;
    ++m;
  }
  catch(std::logic_error& e)
  {
    // do nothing
  }
  for(; m < modsInView.size(); ++m)
  {
    modsInView[m] = {0.0, 0.0};
  }

  // Apply graph
  Jitter.apply(RView, modsInView);

  // Check
  double eps = 1e-15;
  for(std::size_t m = 0; m < RView.size(); ++m)
  {
    CHECK(std::abs(RView[m] - 0.5/QuICC::Math::PI) <= eps);
  }

}

TEST_CASE("Parallel 3D Loop", "[Parallel3DLoop]")
{
  int rank = 0;
  int ranks = 1;
  #ifdef QUICC_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  #endif

  std::string inputFilename = "./simple-3d-loop.mlir";
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
  std::string dist;
  if (rank == 0 && ranks == 1)
  {
    dist = "Serial";
  }
  else if (ranks == 4)
  {
    dist = "Tubular";
  }
  else
  {
    assert(false && "missing setup");
  }
  std::string path = "_refdata/Framework/LoadSplitter/WLFl/";
  std::string id = "103";
  using namespace QuICC::TestSuite;
  auto setup = readDimsAndMeta(path, dist, id);

  // Layouts
  std::array<std::array<std::string, 2>, 3> layOpt;
  layOpt[0] = {"DCCSC3D", "DCCSC3D"};
  layOpt[1] = {"DCCSC3D", "S1CLCSC3D"};
  layOpt[2] = {"DCCSC3D", "DCCSC3D"};

  // Physical space (Stage::PPP and Stage::MPP)
  constexpr std::uint32_t dim = 3;
  std::array<std::vector<std::uint32_t>, dim> pointersPhys {{{}, setup.metaFT.ptr, {}}};
  std::array<std::vector<std::uint32_t>, dim> indicesPhys {{{}, setup.metaFT.idx, {}}};

  // Modal space (aka JW space, Stage::PMM and Stage::MMM, QuICC Stage0)
  std::array<std::vector<std::uint32_t>, dim> pointersMods {{{}, setup.metaJW.ptr, {}}};
  std::array<std::vector<std::uint32_t>, dim> indicesMods {{{}, setup.metaJW.idx, {}}};

  // Store meta stages to pass to Jitter
  std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
  meta.push_back({setup.metaFT.ptr.data(), setup.metaFT.ptr.size()});
  meta.push_back({setup.metaFT.idx.data(), setup.metaFT.idx.size()});
  meta.push_back({setup.metaAL.ptr.data(), setup.metaAL.ptr.size()});
  meta.push_back({setup.metaAL.idx.data(), setup.metaAL.idx.size()});
  meta.push_back({setup.metaJW.ptr.data(), setup.metaJW.ptr.size()});
  meta.push_back({setup.metaJW.idx.data(), setup.metaJW.idx.size()});

  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  using namespace QuICC::Graph;
  Jit<dim> Jitter(std::move(sourceMgr), mem, setup.physDims, setup.modsDims, layOpt, Stage::MMM, Stage::MMM, meta);

  // host mem block
  QuICC::Memory::MemBlock<std::complex<double>> modsIn(setup.modsDims[0]*indicesMods[1].size(), mem.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOut(setup.modsDims[0]*indicesMods[1].size(), mem.get());

  // host view
  C_DCCSC3D_t modsInView({modsIn.data(), modsIn.size()}, {setup.modsDims[0], setup.modsDims[2], setup.modsDims[1]}, pointersMods, indicesMods);
  C_DCCSC3D_t modsOutView({modsOut.data(), modsOut.size()}, {setup.modsDims[0], setup.modsDims[2], setup.modsDims[1]}, pointersMods, indicesMods);

  // set input modes
  std::complex<double> val = {1.0, 0.0};

  for(std::size_t m = 0; m < modsInView.size(); ++m)
  {
    modsInView[m] = {0.0, 0.0};
  }
  try
  {
    // we don't know which rank has mode zero and one
    modsInView(0,0,0);
    // if we didn't throw this rank has mode zero
    modsInView(0,0,0) = val;
    modsInView(1,0,0) = val;
  }
  catch(std::logic_error& e)
  {
    // do nothing
  }

  // Apply graph
  Jitter.apply(modsOutView, modsInView);

  // Check
  double eps = 1e-15;
  for(std::size_t m = 0; m < modsOutView.size(); ++m)
  {
    CHECK(std::abs(modsOutView[m] - modsInView[m]) <= eps);
  }

}

TEST_CASE("Parallel Multi Var 3D Fwd", "[ParallelMultiVar3DFwd]")
{
  int rank = 0;
  int ranks = 1;
  #ifdef QUICC_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  #endif

  std::string inputFilename = "./complex-3d-fwd.mlir";
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
  std::string dist;
  if (rank == 0 && ranks == 1)
  {
    dist = "Serial";
  }
  else if (ranks == 4)
  {
    dist = "Tubular";
  }
  else
  {
    assert(false && "missing setup");
  }
  std::string path = "_refdata/Framework/LoadSplitter/WLFl/";
  std::string id = "103";
  using namespace QuICC::TestSuite;
  auto setup = readDimsAndMeta(path, dist, id);

  // Layouts
  std::array<std::array<std::string, 2>, 3> layOpt;
  layOpt[0] = {"DCCSC3D", "DCCSC3D"};
  layOpt[1] = {"DCCSC3D", "S1CLCSC3D"};
  layOpt[2] = {"DCCSC3D", "DCCSC3D"};

  // Physical space (Stage::PPP and Stage::MPP)
  constexpr std::uint32_t dim = 3;
  std::array<std::vector<std::uint32_t>, dim> pointersPhys {{{}, setup.metaFT.ptr, {}}};
  std::array<std::vector<std::uint32_t>, dim> indicesPhys {{{}, setup.metaFT.idx, {}}};

  // Modal space (aka JW space, Stage::PMM and Stage::MMM, QuICC Stage0)
  std::array<std::vector<std::uint32_t>, dim> pointersMods {{{}, setup.metaJW.ptr, {}}};
  std::array<std::vector<std::uint32_t>, dim> indicesMods {{{}, setup.metaJW.idx, {}}};

  // Store meta stages to pass to Jitter
  std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
  meta.push_back({setup.metaFT.ptr.data(), setup.metaFT.ptr.size()});
  meta.push_back({setup.metaFT.idx.data(), setup.metaFT.idx.size()});
  meta.push_back({setup.metaAL.ptr.data(), setup.metaAL.ptr.size()});
  meta.push_back({setup.metaAL.idx.data(), setup.metaAL.idx.size()});
  meta.push_back({setup.metaJW.ptr.data(), setup.metaJW.ptr.size()});
  meta.push_back({setup.metaJW.idx.data(), setup.metaJW.idx.size()});

  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  using namespace QuICC::Graph;
  Jit<dim> Jitter(std::move(sourceMgr), mem, setup.physDims, setup.modsDims, layOpt, Stage::MMM, Stage::PPP, meta);

  // host mem block
  std::size_t physS = setup.physDims[2]*indicesPhys[1].size();
  QuICC::Memory::MemBlock<double> R(physS, mem.get());
  QuICC::Memory::MemBlock<double> Phi(physS, mem.get());
  QuICC::Memory::MemBlock<double> Theta(physS, mem.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOut(setup.modsDims[0]*indicesMods[1].size(), mem.get());

  // host view
  std::array<std::uint32_t, dim> physDims{setup.physDims[2], setup.physDims[1], setup.physDims[0]};
  R_DCCSC3D_t RView({R.data(), R.size()}, physDims, pointersPhys, indicesPhys);
  R_DCCSC3D_t PhiView({Phi.data(), Phi.size()}, physDims, pointersPhys, indicesPhys);
  R_DCCSC3D_t ThetaView({Theta.data(), Theta.size()}, physDims, pointersPhys, indicesPhys);
  C_DCCSC3D_t modsOutView({modsOut.data(), modsOut.size()}, {setup.modsDims[0], setup.modsDims[2], setup.modsDims[1]}, pointersMods, indicesMods);

  // set input phys
  double val = 1.0;
  for(std::size_t m = 0; m < RView.size(); ++m)
  {
    RView[m] = val;
    PhiView[m] = val;
    ThetaView[m] = val;
  }

  // Apply graph
  Jitter.apply(modsOutView, RView, PhiView, ThetaView);

  // Check
  double eps = 1e-14;
  std::size_t m = 0;
  try
  {
    // we don't know which rank has mode zero
    modsOutView(0,0,0);
    // if we didn't throw this rank has mode zero
    CHECK(std::abs(modsOutView[m].real() - sqrt(2.0)*QuICC::Math::PI) <= eps);
    CHECK(std::abs(modsOutView[m].imag()) <= eps);
    ++m;
  }
  catch(std::logic_error& e)
  {
    // do nothing
  }

  for(; m < modsOutView.size(); ++m)
  {
    CHECK(std::abs(modsOutView[m]) <= eps);
  }
}
