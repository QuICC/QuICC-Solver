
#define CATCH_CONFIG_RUNNER

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


TEST_CASE("One Dimensional Loop Gpu", "[OneDimLoopGpu]")
{
  // Test Graph
  std::string modStr = R"mlir(
    !type_umod = !quiccir.view<5x6x3xf64, "C_DCCSC3D_t">
    !type_uval = !quiccir.view<5x11x3xf64, "R_DCCSC3D_t">
    !type_tumod = tensor<5x6x3xf64, "C_DCCSC3D_t">
    !type_tuval = tensor<5x11x3xf64, "R_DCCSC3D_t">
    func.func @entry(%thisArr: !llvm.ptr<array<2 x ptr>> {llvm.noalias}, %uout: !type_umod, %umod: !type_umod) {
      %tumod = builtin.unrealized_conversion_cast %umod : !type_umod to !type_tumod
      %tuval = quiccir.fr.prj %tumod : !type_tumod -> !type_tuval attributes{implptr = 0 :i64, backend = "cpu"}
      %ret = quiccir.fr.int %tuval : !type_tuval -> !type_tumod attributes{implptr = 1 :i64}
      quiccir.materialize %ret in %uout : (!type_tumod, !type_umod)
      return
    }
  )mlir";

  // Grid dimensions
  constexpr std::uint32_t rank = 3u;
  std::array<std::uint32_t, rank> physDims{11, 3, 5};
  std::array<std::uint32_t, rank> modsDims{6, 3, 5};

  auto memDev = std::make_shared<QuICC::Memory::Cuda::Malloc>();
  QuICC::Graph::Jit<rank> Jitter(modStr, memDev, physDims, modsDims);

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
  QuICC::Profiler::RegionStart<0>("apply-OneDimFourierLoopDealiasGpu");
  Jitter.apply(modsOutViewDev, modsInViewDev);
  QuICC::Profiler::RegionStop<0>("apply-OneDimFourierLoopDealiasGpu");

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
