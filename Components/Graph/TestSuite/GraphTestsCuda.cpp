#include <catch2/catch.hpp>
#include <memory>

// QuICC
#include "Graph/Shims/MlirShims.hpp"
#include "Graph/OpsMap.hpp"
#include "Graph/Jit.hpp"
#include "Graph/Tags.hpp"
#include "Memory/Cpu/NewDelete.hpp"
#include "Memory/Cuda/Malloc.hpp"
#include "Memory/Memory.hpp"
#include "Types/Math.hpp"
#include "Profiler/Interface.hpp"

TEST_CASE("One Dimensional Loop Fourier Gpu", "[OneDimLoopFourierGpu]")
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
  // v012
  std::array<std::uint32_t, rank> physDims{5, 3, 11};
  std::array<std::uint32_t, rank> modsDims{5, 3, 6};

  // View Types
  std::array<std::array<std::string, 2>, 3> layOpt;
  layOpt[0] = {"DCCSC3D", "DCCSC3D"};
  // layOpt[1] = {"DCCSC3DJIK", "S1CLCSC3DJIK"};
  // layOpt[2] = {"DCCSC3DJIK", "DCCSC3DJIK"};

  auto memDev = std::make_shared<QuICC::Memory::Cuda::Malloc>();
  using namespace QuICC::Graph;
  std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
  Jit<rank> Jitter(modStr, memDev, physDims, modsDims, layOpt, Stage::MPP, Stage::MPP, meta);

  // setup metadata
  auto modsM = modsDims[2];
  auto N = physDims[1];
  auto K = physDims[0];
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
    func.func @entry(%tumod: tensor<?x?x?xcomplex<f64>>) -> (tensor<?x?x?xcomplex<f64>>) {
      %tuval = quiccir.al.prj %tumod : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 0 :i64}
      %ret = quiccir.al.int %tuval : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 1 :i64}
      return %ret : tensor<?x?x?xcomplex<f64>>
    }
  )mlir";

  // Grid dimensions
  constexpr std::uint32_t rank = 3u;
  // v012
  std::array<std::uint32_t, rank> physDims{1, 20, 4};
  std::array<std::uint32_t, rank> modsDims{1, 10, 4};

  // View Types
  std::array<std::array<std::string, 2>, 3> layOpt;
  // layOpt[0] = {"DCCSC3D", "DCCSC3D"};
  layOpt[1] = {"DCCSC3DJIK", "S1CLCSC3DJIK"};
  // layOpt[2] = {"DCCSC3DJIK", "DCCSC3DJIK"};

  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  auto memDev = std::make_shared<QuICC::Memory::Cuda::Malloc>();
  using namespace QuICC::Graph;
  std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
  Jit<rank> Jitter(modStr, memDev, physDims, modsDims, layOpt, Stage::MPM, Stage::MPM, meta);

  // setup metadata
  auto M = physDims[2];
  // auto N = physDims[1];
  auto K = physDims[0];
  auto modsM = modsDims[2];
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
    func.func @entry(%tumod: tensor<?x?x?xcomplex<f64>>) -> (tensor<?x?x?xcomplex<f64>>) {
      %tuval = quiccir.jw.prj %tumod : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 0 :i64}
      %ret = quiccir.jw.int %tuval : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 1 :i64}
      return %ret : tensor<?x?x?xcomplex<f64>>
    }
  )mlir";

  // Grid dimensions
  constexpr std::uint32_t rank = 3u;
  // v012
  std::array<std::uint32_t, rank> physDims{3, 4, 1};
  std::array<std::uint32_t, rank> modsDims{2, 4, 1};

  // View Types
  std::array<std::array<std::string, 2>, 3> layOpt;
  // layOpt[0] = {"DCCSC3D", "DCCSC3D"};
  // layOpt[1] = {"DCCSC3DJIK", "S1CLCSC3DJIK"};
  layOpt[2] = {"DCCSC3DJIK", "DCCSC3DJIK"};

  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  auto memDev = std::make_shared<QuICC::Memory::Cuda::Malloc>();
  using namespace QuICC::Graph;
  std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
  Jit<rank> Jitter(modStr, memDev, physDims, modsDims, layOpt, Stage::MMM, Stage::MMM, meta);

  // setup metadata
  auto M = physDims[2];
  auto N = physDims[1];
  // auto K = physDims[0];
  auto modsM = modsDims[2];
  auto modsN = modsDims[1];
  auto modsK = modsDims[0];
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

TEST_CASE("Serial 3D Fwd Gpu", "[Serial3DFwdGpu]")
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
  // v012
  std::array<std::uint32_t, rank> physDims{3, 6, 10};
  // v012
  std::array<std::uint32_t, rank> modsDims{2, 6, 6};

  std::array<std::array<std::string, 2>, 3> layOpt;
  layOpt[0] = {"DCCSC3D", "DCCSC3D"};
  layOpt[1] = {"DCCSC3DJIK", "S1CLCSC3DJIK"};
  layOpt[2] = {"DCCSC3DJIK", "DCCSC3DJIK"};

  // setup metadata
  auto M = physDims[2];
  auto N = physDims[1];
  auto K = physDims[0];
  auto modsM = modsDims[2];
  auto modsN = modsDims[1];
  auto modsK = modsDims[0];

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
  // AL space (Stage::PPM and Stage::MPM)
  using namespace QuICC::Graph;
  auto metaAL = denseTransposePtrAndIdx<C_S1CLCSC3D_t, C_DCCSC3D_t>({modsN, K, modsM});

  // Populate meta for fully populated tensor
  // Modal space (Stage::PMM and Stage::MMM)
  auto metaJW = denseTransposePtrAndIdx<C_DCCSC3D_t, C_S1CLCSC3D_t>({modsK, modsM, modsN});
  std::array<std::vector<std::uint32_t>, rank> pointersMods {{{}, metaJW.ptr, {}}};
  std::array<std::vector<std::uint32_t>, rank> indicesMods {{{}, metaJW.idx, {}}};
  // device meta
  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  auto memDev = std::make_shared<QuICC::Memory::Cuda::Malloc>();
  QuICC::Memory::MemBlock<std::uint32_t> metaFTptrDev(
    ptrPhys.size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> metaFTidxDev(
    idxPhys.size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> metaALptrDev(
    metaAL.ptr.size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> metaALidxDev(
    metaAL.idx.size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> metaJWptrDev(
    metaJW.ptr.size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> metaJWidxDev(
    metaJW.idx.size(), memDev.get());

  // cpu -> gpu index/pointers
  cudaErrChk(cudaMemcpy(metaFTptrDev.data(), ptrPhys.data(),
    ptrPhys.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(metaFTidxDev.data(), idxPhys.data(),
    idxPhys.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(metaALptrDev.data(), metaAL.ptr.data(),
    metaAL.ptr.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(metaALidxDev.data(), metaAL.idx.data(),
    metaAL.idx.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(metaJWptrDev.data(), metaJW.ptr.data(),
    metaJW.ptr.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(metaJWidxDev.data(), metaJW.idx.data(),
    metaJW.idx.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));

  // Store meta stages to pass to Jitter
  std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
  meta.push_back({metaFTptrDev.data(), metaFTptrDev.size()});
  meta.push_back({metaFTidxDev.data(), metaFTidxDev.size()});
  meta.push_back({metaALptrDev.data(), metaALptrDev.size()});
  meta.push_back({metaALidxDev.data(), metaALidxDev.size()});
  meta.push_back({metaJWptrDev.data(), metaJWptrDev.size()});
  meta.push_back({metaJWidxDev.data(), metaJWidxDev.size()});

  // Setup Jitter
  Jit<rank> Jitter(std::move(sourceMgr), memDev, physDims, modsDims, layOpt, Stage::MMM, Stage::PPP, meta);

  // host mem block
  std::size_t physS = M*indicesPhys[1].size();
  QuICC::Memory::MemBlock<double> R(physS, mem.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOut(modsK*metaJW.idx.size(), mem.get());

  // host view
  R_DCCSC3D_t RView({R.data(), R.size()}, {M, N, K}, pointersPhys, indicesPhys);
  C_DCCSC3DJIK_t modsOutView({modsOut.data(), modsOut.size()}, {modsK, modsM, modsN}, pointersMods, indicesMods);

  // set input phys
  double val = 1.0;
  for(std::size_t m = 0; m < RView.size(); ++m)
  {
    RView[m] = val;
  }

  // device block data
  QuICC::Memory::MemBlock<double> RDev(physS, memDev.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOutDev(modsK*metaJW.idx.size(), memDev.get());

  // set device pointers and indices device
  using namespace QuICC::View;
  ViewBase<std::uint32_t> pointersPhysDev[rank];
  pointersPhysDev[1] = ViewBase<std::uint32_t>(metaFTptrDev.data(),
    metaFTptrDev.size());
  ViewBase<std::uint32_t> indicesPhysDev[rank];
  indicesPhysDev[1] = ViewBase<std::uint32_t>(metaFTidxDev.data(),
    metaFTidxDev.size());
  ViewBase<std::uint32_t> pointersModsDev[rank];
  pointersModsDev[1] = ViewBase<std::uint32_t>(metaJWptrDev.data(),
    metaJWptrDev.size());
  ViewBase<std::uint32_t> indicesModsDev[rank];
  indicesModsDev[1] = ViewBase<std::uint32_t>(metaJWidxDev.data(),
    metaJWidxDev.size());

  // device view
  std::array<std::uint32_t, rank> dimPhys = {M, N, K};
  R_DCCSC3D_t RViewDev(RDev.data(), RDev.size(), dimPhys.data(), pointersPhysDev, indicesPhysDev);
  std::array<std::uint32_t, rank> dimMods = {modsK, modsM, modsN};
  C_DCCSC3DJIK_t modsOutViewDev(modsOutDev.data(), modsOutDev.size(), dimMods.data(), pointersModsDev, indicesModsDev);

  // cpu -> gpu data
  cudaErrChk(cudaMemcpy(RDev.data(), R.data(),
    R.size() * sizeof(double), cudaMemcpyHostToDevice));

  // Apply graph
  Jitter.apply(modsOutViewDev, RViewDev);

  // gpu -> cpu
  cudaErrChk(cudaDeviceSynchronize());
  cudaErrChk(cudaMemcpy(modsOut.data(), modsOutDev.data(),
    modsOutDev.size() * sizeof(std::complex<double>),
    cudaMemcpyDeviceToHost));

  // Check
  double eps = 1e-15;
  CHECK(std::abs(modsOutView[0].real() - sqrt(2.0)*QuICC::Math::PI) <= eps);
  CHECK(std::abs(modsOutView[0].imag()) <= eps);
  for(std::size_t m = 1; m < modsOutView.size(); ++m)
  {
    CHECK(std::abs(modsOutView[m]) <= eps);
  }

}

TEST_CASE("Serial 3D Loop Gpu", "[Serial3DLoopGpu]")
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
  // v012
  std::array<std::uint32_t, rank> physDims{6, 10, 10};
  std::array<std::uint32_t, rank> modsDims{3, 7, 6};

  std::array<std::array<std::string, 2>, 3> layOpt;
  layOpt[0] = {"DCCSC3D", "DCCSC3D"};
  layOpt[1] = {"DCCSC3DJIK", "S1CLCSC3DJIK"};
  layOpt[2] = {"DCCSC3DJIK", "DCCSC3DJIK"};

  // setup metadata
  auto N = physDims[1];
  auto K = physDims[0];
  auto modsM = modsDims[2];
  auto modsN = modsDims[1];
  auto modsK = modsDims[0];

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
  // AL space (Stage::PPM and Stage::MPM)
  using namespace QuICC::Graph;
  auto metaAL = denseTransposePtrAndIdx<C_S1CLCSC3DJIK_t, C_DCCSC3DJIK_t>({modsN, K, modsM});

  // Populate meta for fully populated tensor
  // Modal space
  using namespace QuICC::Graph;
  auto metaJW = denseTransposePtrAndIdx<C_DCCSC3DJIK_t, C_S1CLCSC3DJIK_t>({modsK, modsM, modsN});

  // device meta
  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  auto memDev = std::make_shared<QuICC::Memory::Cuda::Malloc>();
  QuICC::Memory::MemBlock<std::uint32_t> metaFTptrDev(
    ptrPhys.size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> metaFTidxDev(
    idxPhys.size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> metaALptrDev(
    metaAL.ptr.size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> metaALidxDev(
    metaAL.idx.size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> metaJWptrDev(
    metaJW.ptr.size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> metaJWidxDev(
    metaJW.idx.size(), memDev.get());

  // cpu -> gpu index/pointers
  cudaErrChk(cudaMemcpy(metaFTptrDev.data(), ptrPhys.data(),
    ptrPhys.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(metaFTidxDev.data(), idxPhys.data(),
    idxPhys.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(metaALptrDev.data(), metaAL.ptr.data(),
    metaAL.ptr.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(metaALidxDev.data(), metaAL.idx.data(),
    metaAL.idx.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(metaJWptrDev.data(), metaJW.ptr.data(),
    metaJW.ptr.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(metaJWidxDev.data(), metaJW.idx.data(),
    metaJW.idx.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));

  // Store meta stages to pass to Jitter
  std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
  meta.push_back({metaFTptrDev.data(), metaFTptrDev.size()});
  meta.push_back({metaFTidxDev.data(), metaFTidxDev.size()});
  meta.push_back({metaALptrDev.data(), metaALptrDev.size()});
  meta.push_back({metaALidxDev.data(), metaALidxDev.size()});
  meta.push_back({metaJWptrDev.data(), metaJWptrDev.size()});
  meta.push_back({metaJWidxDev.data(), metaJWidxDev.size()});

  // Setup Jitter
  using namespace QuICC::Graph;
  Jit<rank> Jitter(std::move(sourceMgr), memDev, physDims, modsDims, layOpt, Stage::MMM, Stage::MMM, meta);

  // host mem block
  std::size_t modsS = modsK * metaJW.idx.size();
  QuICC::Memory::MemBlock<std::complex<double>> modsIn(modsS, mem.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOut(modsS, mem.get());

  // host view
  std::array<std::vector<std::uint32_t>, rank> pointersMods {{{}, metaJW.ptr, {}}};
  std::array<std::vector<std::uint32_t>, rank> indicesMods {{{}, metaJW.idx, {}}};
  C_DCCSC3DJIK_t modsInView({modsIn.data(), modsIn.size()}, {modsK, modsM, modsN}, pointersMods, indicesMods);
  C_DCCSC3DJIK_t modsOutView({modsOut.data(), modsOut.size()}, {modsK, modsM, modsN}, pointersMods, indicesMods);

  // device block data
  QuICC::Memory::MemBlock<std::complex<double>> modsInDev(modsS, memDev.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOutDev(modsS, memDev.get());

  // set device pointers and indices
  using namespace QuICC::View;
  ViewBase<std::uint32_t> pointersDev[rank];
  pointersDev[1] = ViewBase<std::uint32_t>(metaJWptrDev.data(),
    metaJWptrDev.size());
  ViewBase<std::uint32_t> indicesDev[rank];
  indicesDev[1] = ViewBase<std::uint32_t>(metaJWidxDev.data(),
    metaJWidxDev.size());

  // device view
  std::array<std::uint32_t, rank> dimensions = {modsK, modsM, modsN};
  C_DCCSC3DJIK_t modsInViewDev(modsInDev.data(), modsInDev.size(), dimensions.data(), pointersDev, indicesDev);
  C_DCCSC3DJIK_t modsOutViewDev(modsOutDev.data(), modsOutDev.size(), dimensions.data(), pointersDev, indicesDev);

  // set input modes
  std::complex<double> val = {1.0, 0.0};

  for(std::size_t m = 0; m < modsInView.size(); ++m)
  {
    modsInView[m] = {0.0, 0.0};
  }
  modsInView(0,0,0) = val;
  modsInView(1,0,0) = val;

  // cpu -> gpu data
  cudaErrChk(cudaMemcpy(modsInDev.data(), modsIn.data(),
    modsS * sizeof(std::complex<double>), cudaMemcpyHostToDevice));

  // Apply graph
  QuICC::Profiler::RegionStart<0>("apply-SimpleLoopGpu");
  Jitter.apply(modsOutViewDev, modsInViewDev);
  QuICC::Profiler::RegionStop<0>("apply-SimpleLoopGpu");

  // gpu -> cpu
  cudaErrChk(cudaDeviceSynchronize());
  cudaErrChk(cudaMemcpy(modsOut.data(), modsOutDev.data(),
    modsS * sizeof(std::complex<double>),
    cudaMemcpyDeviceToHost));

  // Check
  double eps = 1e-15;
  for(std::size_t m = 0; m < modsOutView.size(); ++m)
  {
    CHECK(std::abs(modsOutView[m] - modsInView[m]) <= eps);
  }

}

TEST_CASE("Serial Multi Var 3D Fwd Gpu", "[SerialMultiVar3DFwdGpu]")
{
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

  // Grid dimensions
  constexpr std::uint32_t rank = 3;
  // v012
  std::array<std::uint32_t, rank> physDims{3, 6, 10};
  std::array<std::uint32_t, rank> modsDims{2, 6, 6};

  std::array<std::array<std::string, 2>, 3> layOpt;
  layOpt[0] = {"DCCSC3D", "DCCSC3D"};
  layOpt[1] = {"DCCSC3DJIK", "S1CLCSC3DJIK"};
  layOpt[2] = {"DCCSC3DJIK", "DCCSC3DJIK"};

  // setup metadata
  auto M = physDims[2];
  auto N = physDims[1];
  auto K = physDims[0];
  auto modsM = modsDims[2];
  auto modsN = modsDims[1];
  auto modsK = modsDims[0];

  std::array<std::uint32_t, 3> inDims {M, N, K};
  std::array<std::uint32_t, 3> outDims {modsK, modsM, modsN};

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
  // AL space (Stage::PPM and Stage::MPM)
  using namespace QuICC::Graph;
  auto metaAL = denseTransposePtrAndIdx<C_S1CLCSC3DJIK_t, C_DCCSC3DJIK_t>({modsN, K, modsM});

  // Populate meta for fully populated tensor
  // Modal space
  using namespace QuICC::Graph;
  auto metaJW = denseTransposePtrAndIdx<C_DCCSC3DJIK_t, C_S1CLCSC3DJIK_t>(outDims);

  // device meta
  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  auto memDev = std::make_shared<QuICC::Memory::Cuda::Malloc>();
  QuICC::Memory::MemBlock<std::uint32_t> metaFTptrDev(
    ptrPhys.size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> metaFTidxDev(
    idxPhys.size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> metaALptrDev(
    metaAL.ptr.size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> metaALidxDev(
    metaAL.idx.size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> metaJWptrDev(
    metaJW.ptr.size(), memDev.get());
  QuICC::Memory::MemBlock<std::uint32_t> metaJWidxDev(
    metaJW.idx.size(), memDev.get());

  // cpu -> gpu index/pointers
  cudaErrChk(cudaMemcpy(metaFTptrDev.data(), ptrPhys.data(),
    ptrPhys.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(metaFTidxDev.data(), idxPhys.data(),
    idxPhys.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(metaALptrDev.data(), metaAL.ptr.data(),
    metaAL.ptr.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(metaALidxDev.data(), metaAL.idx.data(),
    metaAL.idx.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(metaJWptrDev.data(), metaJW.ptr.data(),
    metaJW.ptr.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(metaJWidxDev.data(), metaJW.idx.data(),
    metaJW.idx.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));

  // Store meta stages to pass to Jitter
  std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
  meta.push_back({metaFTptrDev.data(), metaFTptrDev.size()});
  meta.push_back({metaFTidxDev.data(), metaFTidxDev.size()});
  meta.push_back({metaALptrDev.data(), metaALptrDev.size()});
  meta.push_back({metaALidxDev.data(), metaALidxDev.size()});
  meta.push_back({metaJWptrDev.data(), metaJWptrDev.size()});
  meta.push_back({metaJWidxDev.data(), metaJWidxDev.size()});

  // Setup Jitter
  using namespace QuICC::Graph;
  Jit<rank> Jitter(std::move(sourceMgr), memDev, physDims, modsDims, layOpt, Stage::MMM, Stage::PPP, meta);

  // host mem block
  std::size_t physS = M*indicesPhys[1].size();
  std::size_t modsS = modsK*metaJW.idx.size();
  QuICC::Memory::MemBlock<double> R(physS, mem.get());
  QuICC::Memory::MemBlock<double> Phi(physS, mem.get());
  QuICC::Memory::MemBlock<double> Theta(physS, mem.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOut(modsS, mem.get());

  // host view
  std::array<std::vector<std::uint32_t>, rank> pointersMods {{{}, metaJW.ptr, {}}};
  std::array<std::vector<std::uint32_t>, rank> indicesMods {{{}, metaJW.idx, {}}};
  R_DCCSC3D_t RView({R.data(), R.size()}, inDims, pointersPhys, indicesPhys);
  R_DCCSC3D_t PhiView({Phi.data(), Phi.size()}, inDims, pointersPhys, indicesPhys);
  R_DCCSC3D_t ThetaView({Theta.data(), Theta.size()}, inDims, pointersPhys, indicesPhys);
  C_DCCSC3DJIK_t modsOutView({modsOut.data(), modsOut.size()}, outDims, pointersMods, indicesMods);

  // device block data
  QuICC::Memory::MemBlock<double> RDev(physS, memDev.get());
  QuICC::Memory::MemBlock<double> PhiDev(physS, memDev.get());
  QuICC::Memory::MemBlock<double> ThetaDev(physS, memDev.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOutDev(modsS, memDev.get());

  // set device pointers and indices
  using namespace QuICC::View;
  ViewBase<std::uint32_t> pointersPhysDev[rank];
  pointersPhysDev[1] = ViewBase<std::uint32_t>(metaFTptrDev.data(),
    metaFTptrDev.size());
  ViewBase<std::uint32_t> indicesPhysDev[rank];
  indicesPhysDev[1] = ViewBase<std::uint32_t>(metaFTidxDev.data(),
    metaFTidxDev.size());
  ViewBase<std::uint32_t> pointersModsDev[rank];
  pointersModsDev[1] = ViewBase<std::uint32_t>(metaJWptrDev.data(),
    metaJWptrDev.size());
  ViewBase<std::uint32_t> indicesModsDev[rank];
  indicesModsDev[1] = ViewBase<std::uint32_t>(metaJWidxDev.data(),
    metaJWidxDev.size());

  // device view
  R_DCCSC3D_t RViewDev(RDev.data(), RDev.size(), inDims.data(), pointersPhysDev, indicesPhysDev);
  R_DCCSC3D_t PhiViewDev(PhiDev.data(), PhiDev.size(), inDims.data(), pointersPhysDev, indicesPhysDev);
  R_DCCSC3D_t ThetaViewDev(ThetaDev.data(), ThetaDev.size(), inDims.data(), pointersPhysDev, indicesPhysDev);
  C_DCCSC3DJIK_t modsOutViewDev(modsOutDev.data(), modsOutDev.size(), outDims.data(), pointersModsDev, indicesModsDev);

  // set input phys
  double val = 1.0;
  for(std::size_t m = 0; m < RView.size(); ++m)
  {
    RView[m] = val;
    PhiView[m] = val;
    ThetaView[m] = val;
  }

  // cpu -> gpu data
  cudaErrChk(cudaMemcpy(RDev.data(), R.data(),
    physS * sizeof(double), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(PhiDev.data(), Phi.data(),
    physS * sizeof(double), cudaMemcpyHostToDevice));
  cudaErrChk(cudaMemcpy(ThetaDev.data(), Theta.data(),
    physS * sizeof(double), cudaMemcpyHostToDevice));

  // Apply graph
  Jitter.apply(modsOutViewDev, RViewDev, PhiViewDev, ThetaViewDev);

  // gpu -> cpu
  cudaErrChk(cudaDeviceSynchronize());
  cudaErrChk(cudaMemcpy(modsOut.data(), modsOutDev.data(),
    modsS * sizeof(std::complex<double>),
    cudaMemcpyDeviceToHost));

  // Check
  double eps = 1e-15;
  CHECK(std::abs(modsOutView[0].real() - sqrt(2.0)*QuICC::Math::PI) <= eps);
  CHECK(std::abs(modsOutView[0].imag()) <= eps);
  for(std::size_t m = 1; m < modsOutView.size(); ++m)
  {
    CHECK(std::abs(modsOutView[m]) <= eps);
  }
}
