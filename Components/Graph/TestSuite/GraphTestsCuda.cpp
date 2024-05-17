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


TEST_CASE("One Dimensional Loop Fourier Gpu", "[OneDimLoopFourierGpu]")
{
  // Test Graph
  std::string modStr = R"mlir(
    func.func @entry(%tumod: tensor<?x?x?xf64>) -> (tensor<?x?x?xf64>) {
      %tuval = quiccir.fr.prj %tumod : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 0 :i64}
      %ret = quiccir.fr.int %tuval : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 1 :i64}
      return %ret : tensor<?x?x?xf64>
    }
  )mlir";

  // Grid dimensions
  constexpr std::uint32_t rank = 3u;
  std::array<std::uint32_t, rank> physDims{11, 3, 5};
  std::array<std::uint32_t, rank> modsDims{6, 3, 5};

  // View Types
  std::array<std::array<std::string, 2>, 3> layOpt;
  layOpt[0] = {"R_DCCSC3D_t", "C_DCCSC3D_t"};
  // layOpt[1] = {"C_DCCSC3DJIK_t", "C_S1CLCSC3DJIK_t"};
  // layOpt[2] = {"C_DCCSC3DJIK_t", "C_DCCSC3DJIK_t"};

  auto memDev = std::make_shared<QuICC::Memory::Cuda::Malloc>();
  using namespace QuICC::Graph;
  Jit<rank> Jitter(modStr, memDev, physDims, modsDims, layOpt, Stage::MPP, Stage::MPP);

  // setup metadata
  auto modsM = modsDims[0];
  auto N = physDims[1];
  auto K = physDims[2];
  std::array<std::uint32_t, 3> modsDimensions {modsM, N, K};

  std::array<std::vector<std::uint32_t>, 3> pointers {{{},{0, 2, 2, 2, 3, 4},{}}};
  std::array<std::vector<std::uint32_t>, 3> indices {{{},{0, 1, 0, 0},{}}};

  // host mem block
  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  std::size_t modsS = modsM*indices[1].size();
  QuICC::Memory::MemBlock<std::complex<double>> modsIn(modsS, mem.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOut(modsS, mem.get());

  // host view
  using namespace QuICC::Graph;
  C_DCCSC3D_t modsInView({modsIn.data(), modsIn.size()}, modsDimensions, pointers, indices);
  C_DCCSC3D_t modsOutView({modsOut.data(), modsOut.size()}, modsDimensions, pointers, indices);

  // device block
  QuICC::Memory::MemBlock<std::complex<double>> modsInDev(modsS, memDev.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOutDev(modsS, memDev.get());

  QuICC::Memory::MemBlock<std::uint32_t> memBlockPointersDev(
    pointers[1].size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> memBlockIndicesDev(
    indices[1].size(), memDev.get());

  // set device pointers and indice
  using namespace QuICC::View;
  ViewBase<std::uint32_t> pointersDev[rank];
  pointersDev[1] = ViewBase<std::uint32_t>(memBlockPointersDev.data(),
    memBlockPointersDev.size());
  ViewBase<std::uint32_t> indicesDev[rank];
  indicesDev[1] = ViewBase<std::uint32_t>(memBlockIndicesDev.data(),
    memBlockIndicesDev.size());

  // device view
  C_DCCSC3D_t modsInViewDev(modsInDev.data(), modsInDev.size(), modsDimensions.data(), pointersDev, indicesDev);
  C_DCCSC3D_t modsOutViewDev(modsOutDev.data(), modsOutDev.size(), modsDimensions.data(), pointersDev, indicesDev);

  // cpu -> gpu index/pointers
  cudaErrChk(cudaMemcpy(memBlockPointersDev.data(), pointers[1].data(),
    pointers[1].size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(memBlockIndicesDev.data(), indices[1].data(),
    indices[1].size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));

  // set input modes
  std::complex<double> val = {1.0, 0.0};
  for(std::size_t m = 0; m < modsInView.size(); ++m)
  {
    modsInView[m] = val;
  }

  // cpu -> gpu data
  cudaErrChk(cudaMemcpy(modsInDev.data(), modsIn.data(),
    modsS * sizeof(std::complex<double>), cudaMemcpyHostToDevice));

  // Apply graph
  QuICC::Profiler::RegionStart<0>("apply-OneDimFourierLoopGpu");
  Jitter.apply(modsOutViewDev, modsInViewDev);
  QuICC::Profiler::RegionStop<0>("apply-OneDimFourierLoopGpu");

  // gpu -> cpu
  cudaErrChk(cudaDeviceSynchronize());
  cudaErrChk(cudaMemcpy(modsOut.data(), modsOutDev.data(),
    modsS * sizeof(std::complex<double>),
    cudaMemcpyDeviceToHost));

  // Check
  for(std::size_t m = 0; m < modsOutView.size(); ++m)
  {
    CHECK(modsOutView[m] == val);
  }
}

TEST_CASE("One Dimensional Loop Associated Legendre Gpu", "[OneDimLoopALGpu]")
{
  // Test Graph
  // Same setup as transform loop
  std::string modStr = R"mlir(
    func.func @entry(%tumod: tensor<?x?x?xf64>) -> (tensor<?x?x?xf64>) {
      %tuval = quiccir.al.prj %tumod : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 0 :i64}
      %ret = quiccir.al.int %tuval : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 1 :i64}
      return %ret : tensor<?x?x?xf64>
    }
  )mlir";

  // Grid dimensions
  constexpr std::uint32_t rank = 3u;
  std::array<std::uint32_t, rank> physDims{4, 20, 1};
  std::array<std::uint32_t, rank> modsDims{4, 10, 1};

  // View Types
  std::array<std::array<std::string, 2>, 3> layOpt;
  // layOpt[0] = {"R_DCCSC3D_t", "C_DCCSC3D_t"};
  layOpt[1] = {"C_DCCSC3DJIK_t", "C_S1CLCSC3DJIK_t"};
  // layOpt[2] = {"C_DCCSC3DJIK_t", "C_DCCSC3DJIK_t"};

  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  auto memDev = std::make_shared<QuICC::Memory::Cuda::Malloc>();
  using namespace QuICC::Graph;
  Jit<rank> Jitter(modStr, memDev, physDims, modsDims, layOpt, Stage::MPM, Stage::MPM);

  // setup metadata
  auto M = physDims[0];
  // auto N = physDims[1];
  auto K = physDims[2];
  auto modsM = modsDims[0];
  auto modsN = modsDims[1];
  // auto modsK = modsDims[2];
  std::array<std::uint32_t, 3> dimensions {modsN, K, modsM};

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
  C_S1CLCSC3DJIK_t modsInView({modsIn.data(), modsIn.size()}, dimensions, pointers, indices);
  C_S1CLCSC3DJIK_t modsOutView({modsOut.data(), modsOut.size()}, dimensions, pointers, indices);

  // device block
  QuICC::Memory::MemBlock<std::complex<double>> modsInDev(modsS, memDev.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOutDev(modsS, memDev.get());

  QuICC::Memory::MemBlock<std::uint32_t> memBlockPointersDev(
    pointers[1].size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> memBlockIndicesDev(
    indices[1].size(), memDev.get());

  // set device pointers and indice
  using namespace QuICC::View;
  ViewBase<std::uint32_t> pointersDev[rank];
  pointersDev[1] = ViewBase<std::uint32_t>(memBlockPointersDev.data(),
    memBlockPointersDev.size());
  ViewBase<std::uint32_t> indicesDev[rank];
  indicesDev[1] = ViewBase<std::uint32_t>(memBlockIndicesDev.data(),
    memBlockIndicesDev.size());

  // device view
  C_S1CLCSC3DJIK_t modsInViewDev(modsInDev.data(), modsInDev.size(), dimensions.data(), pointersDev, indicesDev);
  C_S1CLCSC3DJIK_t modsOutViewDev(modsOutDev.data(), modsOutDev.size(), dimensions.data(), pointersDev, indicesDev);

  // cpu -> gpu index/pointers
  cudaErrChk(cudaMemcpy(memBlockPointersDev.data(), pointers[1].data(),
    pointers[1].size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(memBlockIndicesDev.data(), indices[1].data(),
    indices[1].size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));

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

  // cpu -> gpu data
  cudaErrChk(cudaMemcpy(modsInDev.data(), modsIn.data(),
    modsS * sizeof(std::complex<double>), cudaMemcpyHostToDevice));

  // Apply graph
  QuICC::Profiler::RegionStart<0>("apply-OneDimALLoopGpu");
  Jitter.apply(modsOutViewDev, modsInViewDev);
  QuICC::Profiler::RegionStop<0>("apply-OneDimALLoopGpu");

  // gpu -> cpu
  cudaErrChk(cudaDeviceSynchronize());
  cudaErrChk(cudaMemcpy(modsOut.data(), modsOutDev.data(),
    modsS * sizeof(std::complex<double>),
    cudaMemcpyDeviceToHost));


  // Check
  double eps = 1e-15;
  for(std::size_t m = 0; m < modsOutView.size(); ++m)
  {
    CHECK(std::abs((modsOutView[m] - modsInView[m]).real()) <= eps);
    CHECK(std::abs((modsOutView[m] - modsInView[m]).imag()) <= eps);
  }
}

TEST_CASE("One Dimensional Loop Worland Gpu", "[OneDimLoopJWGpu]")
{
  // Test Graph
  // Same setup as transform loop
  std::string modStr = R"mlir(
    func.func @entry(%tumod: tensor<?x?x?xf64>) -> (tensor<?x?x?xf64>) {
      %tuval = quiccir.jw.prj %tumod : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 0 :i64}
      %ret = quiccir.jw.int %tuval : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 1 :i64}
      return %ret : tensor<?x?x?xf64>
    }
  )mlir";

  // Grid dimensions
  constexpr std::uint32_t rank = 3u;
  std::array<std::uint32_t, rank> physDims{1, 4, 3};
  std::array<std::uint32_t, rank> modsDims{1, 4, 2};

  // View Types
  std::array<std::array<std::string, 2>, 3> layOpt;
  // layOpt[0] = {"R_DCCSC3D_t", "C_DCCSC3D_t"};
  // layOpt[1] = {"C_DCCSC3DJIK_t", "C_S1CLCSC3DJIK_t"};
  layOpt[2] = {"C_DCCSC3DJIK_t", "C_DCCSC3DJIK_t"};

  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  auto memDev = std::make_shared<QuICC::Memory::Cuda::Malloc>();
  using namespace QuICC::Graph;
  Jit<rank> Jitter(modStr, memDev, physDims, modsDims, layOpt, Stage::MMM, Stage::MMM);

  // setup metadata
  auto M = physDims[0];
  auto N = physDims[1];
  // auto K = physDims[2];
  auto modsM = modsDims[0];
  auto modsN = modsDims[1];
  auto modsK = modsDims[2];
  std::array<std::uint32_t, 3> dimensions {modsK, modsM, modsN};


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
  C_DCCSC3DJIK_t modsInView({modsIn.data(), modsIn.size()}, dimensions, pointers, indices);
  C_DCCSC3DJIK_t modsOutView({modsOut.data(), modsOut.size()}, dimensions, pointers, indices);

  // device block
  QuICC::Memory::MemBlock<std::complex<double>> modsInDev(modsS, memDev.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOutDev(modsS, memDev.get());

  QuICC::Memory::MemBlock<std::uint32_t> memBlockPointersDev(
    pointers[1].size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> memBlockIndicesDev(
    indices[1].size(), memDev.get());

  // set device pointers and indice
  using namespace QuICC::View;
  ViewBase<std::uint32_t> pointersDev[rank];
  pointersDev[1] = ViewBase<std::uint32_t>(memBlockPointersDev.data(),
    memBlockPointersDev.size());
  ViewBase<std::uint32_t> indicesDev[rank];
  indicesDev[1] = ViewBase<std::uint32_t>(memBlockIndicesDev.data(),
    memBlockIndicesDev.size());

  // device view
  C_DCCSC3DJIK_t modsInViewDev(modsInDev.data(), modsInDev.size(), dimensions.data(), pointersDev, indicesDev);
  C_DCCSC3DJIK_t modsOutViewDev(modsOutDev.data(), modsOutDev.size(), dimensions.data(), pointersDev, indicesDev);

  // cpu -> gpu index/pointers
  cudaErrChk(cudaMemcpy(memBlockPointersDev.data(), pointers[1].data(),
    pointers[1].size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(memBlockIndicesDev.data(), indices[1].data(),
    indices[1].size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));

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

  // cpu -> gpu data
  cudaErrChk(cudaMemcpy(modsInDev.data(), modsIn.data(),
    modsS * sizeof(std::complex<double>), cudaMemcpyHostToDevice));

  // Apply graph
  QuICC::Profiler::RegionStart<0>("apply-OneDimJWLoopGpu");
  Jitter.apply(modsOutViewDev, modsInViewDev);
  QuICC::Profiler::RegionStop<0>("apply-OneDimJWLoopGpu");

  // gpu -> cpu
  cudaErrChk(cudaDeviceSynchronize());
  cudaErrChk(cudaMemcpy(modsOut.data(), modsOutDev.data(),
    modsS * sizeof(std::complex<double>),
    cudaMemcpyDeviceToHost));

  // Check
  double eps = 1e-15;
  for(std::size_t m = 0; m < modsOutView.size(); ++m)
  {
    CHECK(std::abs((modsOutView[m] - modsInView[m]).real()) <= eps);
    CHECK(std::abs((modsOutView[m] - modsInView[m]).imag()) <= eps);
  }
}

// TEST_CASE("Serial 3D Loop Gpu", "[Serial3DLoopGpu]")
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

//   // Grid dimensions
//   constexpr std::uint32_t rank = 3;
//   //
//   std::array<std::uint32_t, rank> physDims{10, 10, 6};
//   // M L N
//   std::array<std::uint32_t, rank> modsDims{6, 7, 3};

//   std::array<std::array<std::string, 2>, 3> layOpt;
//   layOpt[0] = {"R_DCCSC3D_t", "C_DCCSC3D_t"};
//   layOpt[1] = {"C_DCCSC3DJIK_t", "C_S1CLCSC3DJIK_t"};
//   layOpt[2] = {"C_DCCSC3DJIK_t", "C_DCCSC3DJIK_t"};

//   auto memDev = std::make_shared<QuICC::Memory::Cuda::Malloc>();
//   using namespace QuICC::Graph;
//   Jit<rank> Jitter(std::move(sourceMgr), memDev, physDims, modsDims, layOpt, Stage::MMM, Stage::MMM);

//   // setup metadata
//   auto modsM = modsDims[0];
//   auto modsN = modsDims[1];
//   auto modsK = modsDims[2];

//   // Populate meta for fully populated tensor
//   // Modal space
//   using namespace QuICC::Graph;
//   auto metaMods = denseTransposePtrAndIdx<C_DCCSC3D_t, C_S1CLCSC3D_t>({modsK, modsM, modsN});
//   std::array<std::vector<std::uint32_t>, rank> pointersMods {{{}, metaMods.ptr, {}}};
//   std::array<std::vector<std::uint32_t>, rank> indicesMods {{{}, metaMods.idx, {}}};

//   // host mem block
//   QuICC::Memory::MemBlock<std::complex<double>> modsIn(modsK*metaMods.idx.size(), mem.get());
//   QuICC::Memory::MemBlock<std::complex<double>> modsOut(modsK*metaMods.idx.size(), mem.get());

//   // host view
//   C_DCCSC3DJIK_t modsInView({modsIn.data(), modsIn.size()}, {modsK, modsM, modsN}, pointersMods, indicesMods);
//   C_DCCSC3DJIK_t modsOutView({modsOut.data(), modsOut.size()}, {modsK, modsM, modsN}, pointersMods, indicesMods);

//   // device block
//   QuICC::Memory::MemBlock<std::complex<double>> modsInDev(modsS, memDev.get());
//   QuICC::Memory::MemBlock<std::complex<double>> modsOutDev(modsS, memDev.get());

//   QuICC::Memory::MemBlock<std::uint32_t> memBlockPointersDev(
//     pointers[1].size(), memDev.get());
//   QuICC::Memory::MemBlock<std::uint32_t> memBlockIndicesDev(
//     indices[1].size(), memDev.get());

//   // set device pointers and indice
//   using namespace QuICC::View;
//   ViewBase<std::uint32_t> pointersDev[rank];
//   pointersDev[1] = ViewBase<std::uint32_t>(memBlockPointersDev.data(),
//     memBlockPointersDev.size());
//   ViewBase<std::uint32_t> indicesDev[rank];
//   indicesDev[1] = ViewBase<std::uint32_t>(memBlockIndicesDev.data(),
//     memBlockIndicesDev.size());

//   // device view
//   C_DCCSC3DJIK_t modsInViewDev(modsInDev.data(), modsInDev.size(), modsDimensions.data(), pointersDev, indicesDev);
//   C_DCCSC3DJIK_t modsOutViewDev(modsOutDev.data(), modsOutDev.size(), modsDimensions.data(), pointersDev, indicesDev);

//   // cpu -> gpu index/pointers
//   cudaErrChk(cudaMemcpy(memBlockPointersDev.data(), pointers[1].data(),
//     pointers[1].size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
//   cudaErrChk(cudaMemcpy(memBlockIndicesDev.data(), indices[1].data(),
//     indices[1].size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));

//   // set input modes
//   std::complex<double> val = {1.0, 0.0};

//   for(std::size_t m = 0; m < modsInView.size(); ++m)
//   {
//     modsInView[m] = {0.0, 0.0};
//   }
//   modsInView(0,0,0) = val;
//   modsInView(1,0,0) = val;

//   // cpu -> gpu data
//   cudaErrChk(cudaMemcpy(modsInDev.data(), modsIn.data(),
//     modsS * sizeof(std::complex<double>), cudaMemcpyHostToDevice));

//   // Apply graph
//   QuICC::Profiler::RegionStart<0>("apply-SimpleLoopGpu");
//   Jitter.apply(modsOutViewDev, modsInViewDev);
//   QuICC::Profiler::RegionStop<0>("apply-SimpleLoopGpu");

//   // gpu -> cpu
//   cudaErrChk(cudaDeviceSynchronize());
//   cudaErrChk(cudaMemcpy(modsOut.data(), modsOutDev.data(),
//     modsS * sizeof(std::complex<double>),
//     cudaMemcpyDeviceToHost));

//   // Check
//   double eps = 1e-15;
//   for(std::size_t m = 0; m < modsOutView.size(); ++m)
//   {
//     CHECK(std::abs(modsOutView[m] - modsInView[m]) <= eps);
//   }

// }
