#include <catch2/catch.hpp>

#include "ViewOps/ViewMemoryUtils.hpp"
#include "ViewOps/Slicewise/Op.hpp"
#include "ViewOps/Slicewise/Functors.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Internal/Casts.hpp"

using namespace QuICC::Memory;
using namespace QuICC::View;

#include "QuICC/Polynomial/Quadrature/WorlandRule.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"


TEST_CASE("Radial Grid Gpu", "RadialGridGpu")
{
   // host mem block
   std::shared_ptr<memory_resource> mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
   std::uint32_t M = 10;
   std::uint32_t N = 4;
   std::uint32_t L = 3;

   using view_t = View<double, DCCSC3D>;
   constexpr std::uint32_t rank = 3;
   std::array<std::uint32_t, rank> dimensions{M, N, L};
   std::array<std::vector<std::uint32_t>, rank> pointers = {
      {{}, {0, 1, 1, 2}, {}}};
   std::array<std::vector<std::uint32_t>, rank> indices = {
      {{}, {0, 1}, {}}};

   std::uint32_t size = M*indices[1].size();

   // host blocks
   // in
   MemBlock<double> memBlockIn(size, mem.get());
   // out and ref
   MemBlock<double> memBlockOut(size, mem.get());
   MemBlock<double> memBlockRef(size, mem.get());

   // host views
   view_t in({memBlockIn.data(), memBlockIn.size()}, dimensions, pointers, indices);
   // out and ref
   view_t out({memBlockOut.data(), memBlockOut.size()}, dimensions, pointers, indices);
   view_t ref({memBlockRef.data(), memBlockRef.size()}, dimensions, pointers, indices);

   // device blocks
   std::shared_ptr<memory_resource> memDev = std::make_shared<QuICC::Memory::Cuda::Malloc>();

   MemBlock<double> memBlockInDev(size, memDev.get());
   // out and ref
   MemBlock<double> memBlockOutDev(size, memDev.get());

   // set device pointers and indices
   QuICC::Memory::MemBlock<std::uint32_t> memBlockPtr(pointers[1].size(), memDev.get());
   QuICC::Memory::MemBlock<std::uint32_t> memBlockIdx(indices[1].size(), memDev.get());

   constexpr std::uint32_t dims = 3;
   ViewBase<std::uint32_t> pointersDev[dims];
   ViewBase<std::uint32_t> indicesDev[dims];

   pointersDev[1] = ViewBase<std::uint32_t>(memBlockPtr.data(), memBlockPtr.size());
   indicesDev[1] = ViewBase<std::uint32_t>(memBlockIdx.data(), memBlockIdx.size());


   // device views
   view_t inDev(memBlockInDev.data(), memBlockInDev.size(), dimensions.data(), pointersDev, indicesDev);
   // out and ref
   view_t outDev(memBlockOutDev.data(), memBlockOutDev.size(), dimensions.data(), pointersDev, indicesDev);

   double scaling = 0.75;

   ::QuICC::Internal::Array igrid;
   ::QuICC::Internal::Array iweights;
   ::QuICC::Polynomial::Quadrature::WorlandRule quad;
   quad.computeQuadrature(igrid, iweights, out.dims()[2]);

   // init col 0, lay 0
   for (std::size_t m = 0; m < M; ++m)
   {
      auto mnl = m;
      in[mnl] = 1.0;
      ref[mnl] = scaling * (in[mnl] * QuICC::Internal::cast(igrid[0]));
   }

   // init col 1, lay 2
   for (std::size_t m = 0; m < M; ++m)
   {
      auto mnl = m + M;
      in[mnl] = 1.0;
      ref[mnl] = scaling * (in[mnl] * QuICC::Internal::cast(igrid[2]));
   }

   // Copy to gpu
   cudaErrChk(cudaMemcpy(inDev.data(), in.data(),
      size * sizeof(typename view_t::ScalarType), cudaMemcpyHostToDevice));
   cudaErrChk(cudaMemcpy(pointersDev[1].data(), pointers[1].data(), pointers[1].size() * sizeof(std::uint32_t),
      cudaMemcpyHostToDevice));
   cudaErrChk(cudaMemcpy(indicesDev[1].data(), indices[1].data(), indices[1].size() * sizeof(std::uint32_t),
      cudaMemcpyHostToDevice));

   // const grid mul op
   using namespace QuICC::Slicewise::Cuda;
   using namespace QuICC::Slicewise;
   auto mulGridOp =
      std::make_unique<Op<2, QuICC::Polynomial::Quadrature::WorlandRule, MulRFunctor<double>,
         view_t, view_t>>(MulRFunctor<double>(scaling), memDev);

   mulGridOp->apply(outDev, inDev);

   // Copy from gpu
   cudaErrChk(cudaMemcpy(out.data(), outDev.data(),
      size * sizeof(typename view_t::ScalarType), cudaMemcpyDeviceToHost));

   // check
   for (std::uint64_t i = 0; i < out.size(); ++i)
   {
      CHECK(out[i] == ref[i]);
   }
}

TEST_CASE("Longitudinal Grid Gpu", "LongitudinalGridGpu")
{
   // host mem block
   std::shared_ptr<memory_resource> mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
   std::uint32_t M = 10;
   std::uint32_t N = 4;
   std::uint32_t L = 3;

   using view_t = View<double, DCCSC3D>;
   constexpr std::uint32_t rank = 3;
   std::array<std::uint32_t, rank> dimensions{M, N, L};
   std::array<std::vector<std::uint32_t>, rank> pointers = {
      {{}, {0, 1, 1, 2}, {}}};
   std::array<std::vector<std::uint32_t>, rank> indices = {
      {{}, {0, 1}, {}}};

   std::uint32_t size = M*indices[1].size();

   // host blocks
   // in
   MemBlock<double> memBlockIn(size, mem.get());
   // out and ref
   MemBlock<double> memBlockOut(size, mem.get());
   MemBlock<double> memBlockRef(size, mem.get());

   // host views
   view_t in({memBlockIn.data(), memBlockIn.size()}, dimensions, pointers, indices);
   // out and ref
   view_t out({memBlockOut.data(), memBlockOut.size()}, dimensions, pointers, indices);
   view_t ref({memBlockRef.data(), memBlockRef.size()}, dimensions, pointers, indices);

   // device blocks
   std::shared_ptr<memory_resource> memDev = std::make_shared<QuICC::Memory::Cuda::Malloc>();

   MemBlock<double> memBlockInDev(size, memDev.get());
   // out and ref
   MemBlock<double> memBlockOutDev(size, memDev.get());

   // set device pointers and indices
   QuICC::Memory::MemBlock<std::uint32_t> memBlockPtr(pointers[1].size(), memDev.get());
   QuICC::Memory::MemBlock<std::uint32_t> memBlockIdx(indices[1].size(), memDev.get());

   constexpr std::uint32_t dims = 3;
   ViewBase<std::uint32_t> pointersDev[dims];
   ViewBase<std::uint32_t> indicesDev[dims];

   pointersDev[1] = ViewBase<std::uint32_t>(memBlockPtr.data(), memBlockPtr.size());
   indicesDev[1] = ViewBase<std::uint32_t>(memBlockIdx.data(), memBlockIdx.size());


   // device views
   view_t inDev(memBlockInDev.data(), memBlockInDev.size(), dimensions.data(), pointersDev, indicesDev);
   // out and ref
   view_t outDev(memBlockOutDev.data(), memBlockOutDev.size(), dimensions.data(), pointersDev, indicesDev);

   double scaling = 0.75;

   ::QuICC::Internal::Array igrid;
   ::QuICC::Internal::Array iweights;
   ::QuICC::Polynomial::Quadrature::LegendreRule quad;
   quad.computeQuadrature(igrid, iweights, out.dims()[1]);
   ::QuICC::Internal::Array itheta = igrid.array().acos();

   // init col 0, lay 0
   for (std::size_t m = 0; m < M; ++m)
   {
      auto mnl = m;
      in[mnl] = 1.0;
      ref[mnl] = scaling * (in[mnl] * QuICC::Internal::cast(itheta[0]));
   }

   // init col 1, lay 2
   for (std::size_t m = 0; m < M; ++m)
   {
      auto mnl = m + M;
      in[mnl] = 1.0;
      ref[mnl] = scaling * (in[mnl] * QuICC::Internal::cast(itheta[1]));
   }

   // Copy to gpu
   cudaErrChk(cudaMemcpy(inDev.data(), in.data(),
      size * sizeof(typename view_t::ScalarType), cudaMemcpyHostToDevice));
   cudaErrChk(cudaMemcpy(pointersDev[1].data(), pointers[1].data(), pointers[1].size() * sizeof(std::uint32_t),
      cudaMemcpyHostToDevice));
   cudaErrChk(cudaMemcpy(indicesDev[1].data(), indices[1].data(), indices[1].size() * sizeof(std::uint32_t),
      cudaMemcpyHostToDevice));

   // const grid mul op
   using namespace QuICC::Slicewise::Cuda;
   using namespace QuICC::Slicewise;
   auto mulGridOp =
      std::make_unique<Op<1, QuICC::Polynomial::Quadrature::LegendreRule, MulSinFunctor<double>,
         view_t, view_t>>(MulSinFunctor<double>(scaling), memDev);

   mulGridOp->apply(outDev, inDev);

   // Copy from gpu
   cudaErrChk(cudaMemcpy(out.data(), outDev.data(),
      size * sizeof(typename view_t::ScalarType), cudaMemcpyDeviceToHost));

   // check
   for (std::uint64_t i = 0; i < out.size(); ++i)
   {
      CHECK(out[i] == ref[i]);
   }
}