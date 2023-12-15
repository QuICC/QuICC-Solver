#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "ViewOps/Reduction/Cpu/Reduction.hpp"
#include "ViewOps/Reduction/Cuda/Reduction.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"

using namespace QuICC::Memory;
using namespace QuICC::View;

TEST_CASE("Reduction Cpu", "ReductionCpu")
{
   using namespace QuICC::Memory;

   // Full data and type
   constexpr size_t M = 4;
   constexpr size_t N = 2;
   constexpr size_t K = 4;

   constexpr size_t S = M * N * 2;
   std::array<double, S> data = {/*k0*/ 1, 2, 3, 4,
      /*k1*/ 5, 6, 7, 8,
      /*k3*/ 9, 10, 11, 12,
      /*k3*/ 13, 14, 15, 16};


   // Akin distributed Fourier op
   // 2D CSC column major matrix (N,K) which elements are 1D dense vectors,
   // i.e a 3D tensor (M,N,K) with fully populated columns
   constexpr std::uint32_t rank = 3;
   std::array<std::uint32_t, rank> dimensions{M, N, K};
   std::array<std::vector<std::uint32_t>, rank> pointers = {
      {{}, {0, 1, 2, 2, 4}, {}}};
   std::array<std::vector<std::uint32_t>, rank> indices = {
      {{}, {0, 0, 0, 1}, {}}};
   View<double, DCCSC3D> view3D(data, dimensions, pointers, indices);

   // Reduced data and ref
   constexpr size_t SR = N * 2;
   std::array<double, SR> data2D;
   std::array<double, SR> data2Dref = {10, 26, 42, 58};

   constexpr std::uint32_t rank2D = 2;
   std::array<std::uint32_t, rank2D> dimensions2D{N, K};
   std::array<std::vector<std::uint32_t>, rank2D> pointers2D = {
      {pointers[1], {}}};
   std::array<std::vector<std::uint32_t>, rank2D> indices2D = {
      {indices[1], {}}};

   View<double, CSC> view2D(data2D, dimensions2D, pointers2D, indices2D);

   using namespace QuICC::Reduction::Cpu;
   auto reductorOp =
      std::make_unique<Op<View<double, CSC>, View<double, DCCSC3D>, 0u>>();

   reductorOp->apply(view2D, view3D);

   // check
   for (std::uint64_t i = 0; i < data2D.size(); ++i)
   {
      CHECK(data2D[i] == data2Dref[i]);
   }
}

#ifdef QUICC_HAS_CUDA_BACKEND

TEST_CASE("Reduction Gpu", "ReductionGpu")
{
   using namespace QuICC::Memory;

   // Full data and type
   constexpr size_t M = 4;
   constexpr size_t N = 2;
   constexpr size_t K = 4;

   constexpr size_t S = M * N * 2;
   std::array<double, S> data = {/*k0*/ 1, 2, 3, 4,
      /*k1*/ 5, 6, 7, 8,
      /*k3*/ 9, 13,
      /*k3*/ 10, 14,
      /*k3*/ 11, 15,
      /*k3*/ 12, 16};


   // Akin distributed Fourier op
   // 2D CSC column major matrix (N,K) which elements are 1D dense vectors,
   // i.e a 3D tensor (M,N,K) with fully populated columns
   constexpr std::uint32_t rank = 3;
   std::array<std::uint32_t, rank> dimensions{M, N, K};
   std::array<std::vector<std::uint32_t>, rank> pointers = {
      {{}, {0, 1, 2, 2, 4}, {}}};
   std::array<std::vector<std::uint32_t>, rank> indices = {
      {{}, {0, 0, 0, 1}, {}}};
   using view3D_t = View<double, DCCSC3DJIK>;
   view3D_t view3D(data, dimensions, pointers, indices);

   // Reduced data and ref
   constexpr size_t SR = N * 2;
   std::array<double, SR> data2D;
   std::array<double, SR> data2Dref = {10, 26, 42, 58};

   constexpr std::uint32_t rank2D = 2;
   std::array<std::uint32_t, rank2D> dimensions2D{N, K};
   std::array<std::vector<std::uint32_t>, rank2D> pointers2D = {
      {pointers[1], {}}};
   std::array<std::vector<std::uint32_t>, rank2D> indices2D = {
      {indices[1], {}}};
   using view2D_t = View<double, CSC>;
   view2D_t view2D(data2D, dimensions2D, pointers2D, indices2D);


   // device mem 3D blocks
   auto memDev = std::make_shared<QuICC::Memory::Cuda::Malloc>();
   QuICC::Memory::MemBlock<double> memBlock3DDev(S, memDev.get());
   QuICC::Memory::MemBlock<std::uint32_t> memBlock3DDevPointers(
      pointers[1].size(), memDev.get());
   QuICC::Memory::MemBlock<std::uint32_t> memBlock3DDevIndices(
      indices[1].size(), memDev.get());

   // set device pointers and indices 3D
   ViewBase<std::uint32_t> pointersDev[rank];
   pointersDev[1] = ViewBase<std::uint32_t>(memBlock3DDevPointers.data(),
      memBlock3DDevPointers.size());
   ViewBase<std::uint32_t> indicesDev[rank];
   indicesDev[1] = ViewBase<std::uint32_t>(memBlock3DDevIndices.data(),
      memBlock3DDevIndices.size());

   // set device views 3D
   view3D_t view3DDev(memBlock3DDev.data(), memBlock3DDev.size(),
      dimensions.data(), pointersDev, indicesDev);

   // device mem 2D blocks
   QuICC::Memory::MemBlock<double> memBlock2DDev(view2D.size(), memDev.get());
   QuICC::Memory::MemBlock<std::uint32_t> memBlock2DDevPointers(
      pointers2D[0].size(), memDev.get());
   QuICC::Memory::MemBlock<std::uint32_t> memBlock2DDevIndices(
      indices2D[0].size(), memDev.get());

   // set device pointers and indices 2D
   ViewBase<std::uint32_t> pointers2DDev[rank2D];
   pointers2DDev[0] = ViewBase<std::uint32_t>(memBlock2DDevPointers.data(),
      memBlock2DDevPointers.size());
   ViewBase<std::uint32_t> indices2DDev[rank2D];
   indices2DDev[0] = ViewBase<std::uint32_t>(memBlock2DDevIndices.data(),
      memBlock2DDevIndices.size());

   // set device views 2D
   view2D_t view2DDev(memBlock2DDev.data(), memBlock2DDev.size(),
      dimensions2D.data(), pointers2DDev, indices2DDev);

   // cpu -> gpu 3D
   cudaErrChk(cudaMemcpy(view3DDev.data(), view3D.data(),
      S * sizeof(typename view3D_t::ScalarType), cudaMemcpyHostToDevice));
   cudaErrChk(cudaMemcpy(view3DDev.pointers()[1].data(), pointers[1].data(),
      pointers[1].size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
   cudaErrChk(cudaMemcpy(view3DDev.indices()[1].data(), indices[1].data(),
      indices[1].size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));

   // cpu -> gpu 2D
   cudaErrChk(cudaMemcpy(view2DDev.data(), view2D.data(),
      view2D.size() * sizeof(typename view2D_t::ScalarType),
      cudaMemcpyHostToDevice));
   cudaErrChk(cudaMemcpy(view2DDev.pointers()[0].data(), pointers2D[0].data(),
      pointers2D[0].size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));
   cudaErrChk(cudaMemcpy(view2DDev.indices()[0].data(), indices2D[0].data(),
      indices2D[0].size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice));


   using namespace QuICC::Reduction::Cuda;
   auto reductorOp =
      std::make_unique<Op<View<double, CSC>, View<double, DCCSC3DJIK>, 0u>>(
         memDev);

   reductorOp->apply(view2DDev, view3DDev);

   // gpu -> cpu
   cudaErrChk(cudaMemcpy(view2D.data(), view2DDev.data(),
      view2D.size() * sizeof(typename view2D_t::ScalarType),
      cudaMemcpyDeviceToHost));

   // check
   for (std::uint64_t i = 0; i < data2D.size(); ++i)
   {
      CHECK(data2D[i] == data2Dref[i]);
   }
}


#endif
