#include <catch2/catch.hpp>
#include <memory>

// QuICC
#include "Graph/Shims/MlirShims.hpp"
#include "Graph/BackendsMap.hpp"
#include "Graph/OpsMap.hpp"
#include "Graph/Jit.hpp"
#include "Memory/Cpu/NewDelete.hpp"
#include "Memory/Cuda/Malloc.hpp"
#include "Memory/Memory.hpp"

#include "Types/Math.hpp"
#include "Profiler/Interface.hpp"

TEST_CASE("One Dimensional Loop Fourier", "[OneDimLoopFourier]")
{
  // Test Graph
  std::string modStr = R"mlir(
    func.func @entry(%tumod: tensor<?x?x?xcomplex<f64>>) -> (tensor<?x?x?xcomplex<f64>>) {
      %tuval = quiccir.fr.prj %tumod : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xf64> attributes{implptr = 0 :i64}
      %ret = quiccir.fr.int %tuval : tensor<?x?x?xf64> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 1 :i64}
      return %ret : tensor<?x?x?xcomplex<f64>>
    }
  )mlir";

  // Grid dimensions
  constexpr std::uint32_t rank = 3u;
  std::array<std::uint32_t, rank> physDims{11, 3, 5};
  std::array<std::uint32_t, rank> modsDims{6, 3, 5};

  // View Types
  std::array<std::array<std::string, 2>, 3> layOpt;
  layOpt[0] = {"DCCSC3D", "DCCSC3D"};
  // layOpt[1] = {"DCCSC3D", "S1CLCSC3D"};
  // layOpt[2] = {"DCCSC3D", "DCCSC3D"};

  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  using namespace QuICC::Graph;
  std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
  Jit<rank> Jitter(modStr, mem, physDims, modsDims, layOpt, Stage::MPP, Stage::MPP, meta);

  // setup metadata
  auto modsM = modsDims[0];
  auto N = physDims[1];
  auto K = physDims[2];
  std::array<std::uint32_t, 3> modsDimensions {modsM, N, K};

  std::array<std::vector<std::uint32_t>, 3> pointers {{{},{0, 2, 2, 2, 3, 4},{}}};
  std::array<std::vector<std::uint32_t>, 3> indices {{{},{0, 1, 0, 0},{}}};

  // host mem block
  std::size_t modsS = modsM*indices[1].size();
  QuICC::Memory::MemBlock<std::complex<double>> modsIn(modsS, mem.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOut(modsS, mem.get());

  // host view
  using namespace QuICC::Graph;
  C_DCCSC3D_t modsInView({modsIn.data(), modsIn.size()}, modsDimensions, pointers, indices);
  C_DCCSC3D_t modsOutView({modsOut.data(), modsOut.size()}, modsDimensions, pointers, indices);

  // set input modes
  std::complex<double> val = {1.0, 0.0};
  for(std::size_t m = 0; m < modsInView.size(); ++m)
  {
    modsInView[m] = val;
  }

  // Apply graph
  Jitter.apply(modsOutView, modsInView);

  // Check
  for(std::size_t m = 0; m < modsOutView.size(); ++m)
  {
    CHECK(modsOutView[m] == val);
  }
}

TEST_CASE("One Dimensional Loop Associated Legendre", "[OneDimLoopAL]")
{
  // Test Graph
  // Same setup as transform loop
  std::string modStr = R"mlir(
    func.func @entry(%tumod: tensor<?x?x?xcomplex<f64>>) -> (tensor<?x?x?xcomplex<f64>>) {
      %tuval = quiccir.al.prj %tumod : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 0 :i64}
      %ret = quiccir.al.int %tuval : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 1 :i64}
      return %ret : tensor<?x?x?xcomplex<f64>>
    }
  )mlir";

  // Grid dimensions
  constexpr std::uint32_t rank = 3u;
  std::array<std::uint32_t, rank> physDims{4, 20, 1};
  std::array<std::uint32_t, rank> modsDims{4, 10, 1};

  // View Types
  std::array<std::array<std::string, 2>, 3> layOpt;
  // layOpt[0] = {"DCCSC3D", "DCCSC3D"};
  layOpt[1] = {"DCCSC3D", "S1CLCSC3D"};
  // layOpt[2] = {"DCCSC3D", "DCCSC3D"};

  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  using namespace QuICC::Graph;
  std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
  Jit<rank> Jitter(modStr, mem, physDims, modsDims, layOpt, Stage::MPM, Stage::MPM, meta);

  // setup metadata
  auto M = physDims[0];
  // auto N = physDims[1];
  auto K = physDims[2];
  auto modsM = modsDims[0];
  auto modsN = modsDims[1];
  // auto modsK = modsDims[2];

  // Populate meta for fully populated tensor

  std::vector<std::uint32_t> ptr(M+1);
  std::vector<std::uint32_t> idx(M*K);
  ptr[0] = 0;
  for (std::size_t i = 1; i < ptr.size(); ++i) {
      ptr[i] = ptr[i-1]+K;
  }
  for (std::size_t i = 0; i < idx.size(); ++i) {
      idx[i] = i % K;
  }
  std::array<std::vector<std::uint32_t>, rank> pointers {{{}, ptr, {}}};
  std::array<std::vector<std::uint32_t>, rank> indices {{{}, idx, {}}};

  // host mem block
  std::size_t modsS = modsN + modsN-1 + modsN-2 + modsN-3;
  QuICC::Memory::MemBlock<std::complex<double>> modsIn(modsS, mem.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOut(modsS, mem.get());

  // host view
  using namespace QuICC::Graph;
  C_S1CLCSC3D_t modsInView({modsIn.data(), modsIn.size()}, {modsN, K, modsM}, pointers, indices);
  C_S1CLCSC3D_t modsOutView({modsOut.data(), modsOut.size()}, {modsN, K, modsM}, pointers, indices);

  // set input modes
  std::complex<double> val = {1.0, -1.0};
  for(std::size_t m = 0; m < modsInView.size(); ++m)
  {
    modsInView[m] = {0.0, 0.0};
  }
  modsInView(0, 0, 0) = val;
  modsInView(0, 0, 1) = val;
  modsInView(0, 0, 2) = val;
  modsInView(0, 0, 3) = val;

  // Apply graph
  Jitter.apply(modsOutView, modsInView);

  // Check
  double eps = 1e-15;
  for(std::size_t m = 0; m < modsOutView.size(); ++m)
  {
    CHECK(std::abs((modsOutView[m] - modsInView[m]).real()) <= eps);
    CHECK(std::abs((modsOutView[m] - modsInView[m]).imag()) <= eps);
  }
}

TEST_CASE("One Dimensional Loop Worland", "[OneDimLoopJW]")
{
  // Test Graph
  // Same setup as transform loop
  std::string modStr = R"mlir(
    func.func @entry(%tumod: tensor<?x?x?xcomplex<f64>>) -> (tensor<?x?x?xcomplex<f64>>) {
      %tuval = quiccir.jw.prj %tumod : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 0 :i64}
      %ret = quiccir.jw.int %tuval : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 1 :i64}
      return %ret : tensor<?x?x?xcomplex<f64>>
    }
  )mlir";

  // Grid dimensions
  constexpr std::uint32_t rank = 3u;
  std::array<std::uint32_t, rank> physDims{1, 4, 3};
  std::array<std::uint32_t, rank> modsDims{1, 4, 2};

  // View Types
  std::array<std::array<std::string, 2>, 3> layOpt;
  // layOpt[0] = {"DCCSC3D", "DCCSC3D"};
  // layOpt[1] = {"DCCSC3D", "S1CLCSC3D"};
  layOpt[2] = {"DCCSC3D", "DCCSC3D"};

  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  using namespace QuICC::Graph;
  std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
  Jit<rank> Jitter(modStr, mem, physDims, modsDims, layOpt, Stage::MMM, Stage::MMM, meta);

  // setup metadata
  auto M = physDims[0];
  auto N = physDims[1];
  // auto K = physDims[2];
  auto modsM = modsDims[0];
  auto modsN = modsDims[1];
  auto modsK = modsDims[2];

  // Populate meta for fully populated tensor
  std::vector<std::uint32_t> ptr(N+1);
  std::vector<std::uint32_t> idx(N*M);
  ptr[0] = 0;
  for (std::size_t i = 1; i < ptr.size(); ++i) {
      ptr[i] = ptr[i-1]+M;
  }
  for (std::size_t i = 0; i < idx.size(); ++i) {
      idx[i] = i % M;
  }
  std::array<std::vector<std::uint32_t>, rank> pointers {{{}, ptr, {}}};
  std::array<std::vector<std::uint32_t>, rank> indices {{{}, idx, {}}};

  // host mem block
  std::size_t modsS = modsK * idx.size();
  QuICC::Memory::MemBlock<std::complex<double>> modsIn(modsS, mem.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOut(modsS, mem.get());

  // host view
  using namespace QuICC::Graph;
  C_DCCSC3D_t modsInView({modsIn.data(), modsIn.size()}, {modsK, modsM, modsN}, pointers, indices);
  C_DCCSC3D_t modsOutView({modsOut.data(), modsOut.size()}, {modsK, modsM, modsN}, pointers, indices);

  // set input modes
  std::complex<double> val = {1.0, 0.0};
  for(std::size_t m = 0; m < modsInView.size(); ++m)
  {
    modsInView[m] = {0.0, 0.0};
  }
  modsInView(0, 0, 0) = val;
  modsInView(0, 0, 1) = val;
  modsInView(0, 0, 2) = val;
  modsInView(0, 0, 3) = val;

  // Apply graph
  Jitter.apply(modsOutView, modsInView);

  // Check
  double eps = 1e-15;
  for(std::size_t m = 0; m < modsOutView.size(); ++m)
  {
    CHECK(std::abs((modsOutView[m] - modsInView[m]).real()) <= eps);
    CHECK(std::abs((modsOutView[m] - modsInView[m]).imag()) <= eps);
  }
}

TEST_CASE("Serial 3D Fwd", "[Serial3DFwd]")
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

  // Grid dimensions
  constexpr std::uint32_t rank = 3;
  std::array<std::uint32_t, rank> physDims{10, 6, 3};
  std::array<std::uint32_t, rank> modsDims{6, 6, 2};

  std::array<std::array<std::string, 2>, 3> layOpt;
  layOpt[0] = {"DCCSC3D", "DCCSC3D"};
  layOpt[1] = {"DCCSC3D", "S1CLCSC3D"};
  layOpt[2] = {"DCCSC3D", "DCCSC3D"};

  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();

  // setup metadata
  auto M = physDims[0];
  auto N = physDims[1];
  auto K = physDims[2];
  auto modsM = modsDims[0];
  auto modsN = modsDims[1];
  auto modsK = modsDims[2];

  // Populate meta for fully populated tensor
  // Physical space (Stage::PPP and Stage::MPP)
  std::vector<std::uint32_t> ptrPhys(K+1);
  std::vector<std::uint32_t> idxPhys(K*N);
  ptrPhys[0] = 0;
  for (std::size_t i = 1; i < ptrPhys.size(); ++i) {
      ptrPhys[i] = ptrPhys[i-1]+N;
  }
  for (std::size_t i = 0; i < idxPhys.size(); ++i) {
      idxPhys[i] = i % N;
  }
  std::array<std::vector<std::uint32_t>, rank> pointersPhys {{{}, ptrPhys, {}}};
  std::array<std::vector<std::uint32_t>, rank> indicesPhys {{{}, idxPhys, {}}};

  // Populate meta for fully populated tensor
  // AL space (Stage::PPM and Stage::MPM
  using namespace QuICC::Graph;
  auto metaAL = denseTransposePtrAndIdx<C_S1CLCSC3D_t, C_DCCSC3D_t>({modsN, K, modsM});

  // Populate meta for fully populated tensor
  // Modal space (Stage::PMM and Stage::MMM)
  auto metaJW = denseTransposePtrAndIdx<C_DCCSC3D_t, C_S1CLCSC3D_t>({modsK, modsM, modsN});
  std::array<std::vector<std::uint32_t>, rank> pointersMods {{{}, metaJW.ptr, {}}};
  std::array<std::vector<std::uint32_t>, rank> indicesMods {{{}, metaJW.idx, {}}};

  // Store meta stages to pass to Jitter
  std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
  meta.push_back({nullptr, 0});
  meta.push_back({nullptr, 0});
  meta.push_back({metaAL.ptr.data(), metaAL.ptr.size()});
  meta.push_back({metaAL.idx.data(), metaAL.idx.size()});
  meta.push_back({metaJW.ptr.data(), metaJW.ptr.size()});
  meta.push_back({metaJW.idx.data(), metaJW.idx.size()});

  // Setup Jitter
  using namespace QuICC::Graph;
  Jit<rank> Jitter(std::move(sourceMgr), mem, physDims, modsDims, layOpt, Stage::MMM, Stage::PPP, meta);

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
  double eps = 1e-15;
  CHECK(std::abs(modsOutView[0].real() - sqrt(2.0)*QuICC::Math::PI) <= eps);
  CHECK(std::abs(modsOutView[0].imag()) <= eps);
  for(std::size_t m = 1; m < modsOutView.size(); ++m)
  {
    CHECK(std::abs(modsOutView[m]) <= eps);
  }

}

TEST_CASE("Serial 3D Bwd", "[Serial3DBwd]")
{
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

  // Grid dimensions
  constexpr std::uint32_t rank = 3;
  std::array<std::uint32_t, rank> physDims{10, 6, 3};
  std::array<std::uint32_t, rank> modsDims{6, 6, 2};

  std::array<std::array<std::string, 2>, 3> layOpt;
  layOpt[0] = {"DCCSC3D", "DCCSC3D"};
  layOpt[1] = {"DCCSC3D", "S1CLCSC3D"};
  layOpt[2] = {"DCCSC3D", "DCCSC3D"};

  // setup metadata
  auto M = physDims[0];
  auto N = physDims[1];
  auto K = physDims[2];
  auto modsM = modsDims[0];
  auto modsN = modsDims[1];
  auto modsK = modsDims[2];

  // Populate meta for fully populated tensor
  // Physical space
  std::vector<std::uint32_t> ptrPhys(K+1);
  std::vector<std::uint32_t> idxPhys(K*N);
  ptrPhys[0] = 0;
  for (std::size_t i = 1; i < ptrPhys.size(); ++i) {
      ptrPhys[i] = ptrPhys[i-1]+N;
  }
  for (std::size_t i = 0; i < idxPhys.size(); ++i) {
      idxPhys[i] = i % N;
  }
  std::array<std::vector<std::uint32_t>, rank> pointersPhys {{{}, ptrPhys, {}}};
  std::array<std::vector<std::uint32_t>, rank> indicesPhys {{{}, idxPhys, {}}};

  // Populate meta for fully populated tensor
  // AL space (Stage::PPM and Stage::MPM
  using namespace QuICC::Graph;
  auto metaAL = denseTransposePtrAndIdx<C_S1CLCSC3D_t, C_DCCSC3D_t>({modsN, K, modsM});

  // Populate meta for fully populated tensor
  // Modal space
  auto metaJW = denseTransposePtrAndIdx<C_DCCSC3D_t, C_S1CLCSC3D_t>({modsK, modsM, modsN});
  std::array<std::vector<std::uint32_t>, rank> pointersMods {{{}, metaJW.ptr, {}}};
  std::array<std::vector<std::uint32_t>, rank> indicesMods {{{}, metaJW.idx, {}}};

  // Store meta stages to pass to Jitter
  std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
  meta.push_back({ptrPhys.data(), ptrPhys.size()});
  meta.push_back({idxPhys.data(), idxPhys.size()});
  meta.push_back({metaAL.ptr.data(), metaAL.ptr.size()});
  meta.push_back({metaAL.idx.data(), metaAL.idx.size()});
  meta.push_back({metaJW.ptr.data(), metaJW.ptr.size()});
  meta.push_back({metaJW.idx.data(), metaJW.idx.size()});

  // Setup Jitter
  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  using namespace QuICC::Graph;
  Jit<rank> Jitter(std::move(sourceMgr), mem, physDims, modsDims, layOpt, Stage::PPP, Stage::MMM, meta);

  // host mem block
  std::size_t physS = M*indicesPhys[1].size();
  QuICC::Memory::MemBlock<double> R(physS, mem.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsIn(modsK*metaJW.idx.size(), mem.get());

  // host view
  R_DCCSC3D_t RView({R.data(), R.size()}, physDims, pointersPhys, indicesPhys);
  C_DCCSC3D_t modsInView({modsIn.data(), modsIn.size()}, {modsK, modsM, modsN}, pointersMods, indicesMods);

  // set input modes
  std::complex<double> val = {sqrt(2.0)/2.0, 0.0};
  modsInView[0] = val;
  for(std::size_t m = 1; m < modsInView.size(); ++m)
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

TEST_CASE("Serial 3D Loop", "[Serial3DLoop]")
{
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

  // Grid dimensions
  constexpr std::uint32_t rank = 3;
  //
  std::array<std::uint32_t, rank> physDims{10, 10, 6};
  // M L N
  std::array<std::uint32_t, rank> modsDims{6, 7, 3};

  std::array<std::array<std::string, 2>, 3> layOpt;
  layOpt[0] = {"DCCSC3D", "DCCSC3D"};
  layOpt[1] = {"DCCSC3D", "S1CLCSC3D"};
  layOpt[2] = {"DCCSC3D", "DCCSC3D"};

  // setup metadata
  auto M = physDims[0];
  auto N = physDims[1];
  auto K = physDims[2];
  auto modsM = modsDims[0];
  auto modsN = modsDims[1];
  auto modsK = modsDims[2];

  // Populate meta for fully populated tensor
  // Physical space
  std::vector<std::uint32_t> ptrPhys(K+1);
  std::vector<std::uint32_t> idxPhys(K*N);
  ptrPhys[0] = 0;
  for (std::size_t i = 1; i < ptrPhys.size(); ++i) {
      ptrPhys[i] = ptrPhys[i-1]+N;
  }
  for (std::size_t i = 0; i < idxPhys.size(); ++i) {
      idxPhys[i] = i % N;
  }
  std::array<std::vector<std::uint32_t>, rank> pointersPhys {{{}, ptrPhys, {}}};
  std::array<std::vector<std::uint32_t>, rank> indicesPhys {{{}, idxPhys, {}}};

  // Populate meta for fully populated tensor
  // AL space (Stage::PPM and Stage::MPM
  using namespace QuICC::Graph;
  auto metaAL = denseTransposePtrAndIdx<C_S1CLCSC3D_t, C_DCCSC3D_t>({modsN, K, modsM});

  // Populate meta for fully populated tensor
  // Modal space
  using namespace QuICC::Graph;
  auto metaJW = denseTransposePtrAndIdx<C_DCCSC3D_t, C_S1CLCSC3D_t>({modsK, modsM, modsN});
  std::array<std::vector<std::uint32_t>, rank> pointersMods {{{}, metaJW.ptr, {}}};
  std::array<std::vector<std::uint32_t>, rank> indicesMods {{{}, metaJW.idx, {}}};

  // Store meta stages to pass to Jitter
  std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
  meta.push_back({ptrPhys.data(), ptrPhys.size()});
  meta.push_back({idxPhys.data(), idxPhys.size()});
  meta.push_back({metaAL.ptr.data(), metaAL.ptr.size()});
  meta.push_back({metaAL.idx.data(), metaAL.idx.size()});
  meta.push_back({metaJW.ptr.data(), metaJW.ptr.size()});
  meta.push_back({metaJW.idx.data(), metaJW.idx.size()});

  // Setup Jitter
  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  using namespace QuICC::Graph;
  Jit<rank> Jitter(std::move(sourceMgr), mem, physDims, modsDims, layOpt, Stage::MMM, Stage::MMM, meta);

  // host mem block
  QuICC::Memory::MemBlock<std::complex<double>> modsIn(modsK*metaJW.idx.size(), mem.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOut(modsK*metaJW.idx.size(), mem.get());

  // host view
  C_DCCSC3D_t modsInView({modsIn.data(), modsIn.size()}, {modsK, modsM, modsN}, pointersMods, indicesMods);
  C_DCCSC3D_t modsOutView({modsOut.data(), modsOut.size()}, {modsK, modsM, modsN}, pointersMods, indicesMods);

  // set input modes
  std::complex<double> val = {1.0, 0.0};

  for(std::size_t m = 0; m < modsInView.size(); ++m)
  {
    modsInView[m] = {0.0, 0.0};
  }
  modsInView(0,0,0) = val;
  modsInView(1,0,0) = val;

  // Apply graph
  Jitter.apply(modsOutView, modsInView);

  // Check
  double eps = 1e-15;
  for(std::size_t m = 0; m < modsOutView.size(); ++m)
  {
    CHECK(std::abs(modsOutView[m] - modsInView[m]) <= eps);
  }

}

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
//   layOpt[0] = {"DCCSC3D", "DCCSC3D"};
//   layOpt[1] = {"DCCSC3D", "S1CLCSC3D"};
//   layOpt[2] = {"DCCSC3D", "DCCSC3D"};

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
