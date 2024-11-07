#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "ViewOps/Slicewise/Op.hpp"
#include "ViewOps/Slicewise/Functors.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"

#include "QuICC/Polynomial/Quadrature/WorlandRule.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"

using namespace QuICC::Memory;
using namespace QuICC::View;

TEST_CASE("Radial Grid", "RadialGrid")
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

   // in
   MemBlock<double> memBlockIn(size, mem.get());
   // out and ref
   MemBlock<double> memBlockOut(size, mem.get());
   MemBlock<double> memBlockRef(size, mem.get());

   // views
   view_t in({memBlockIn.data(), memBlockIn.size()}, dimensions, pointers, indices);
   // out and ref
   view_t out({memBlockOut.data(), memBlockOut.size()}, dimensions, pointers, indices);
   view_t ref({memBlockRef.data(), memBlockRef.size()}, dimensions, pointers, indices);

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

   // const grid mul op
   using namespace QuICC::Slicewise::Cpu;
   using namespace QuICC::Slicewise;
   auto mulGridOp =
      std::make_unique<Op<2, ::QuICC::Polynomial::Quadrature::WorlandRule, MulRFunctor<double>,
         view_t, view_t>>(MulRFunctor<double>(scaling), mem);

   mulGridOp->apply(out, in);

   // check
   for (std::uint64_t i = 0; i < out.size(); ++i)
   {
      CHECK(out[i] == ref[i]);
   }
}

TEST_CASE("Longitudinal Grid", "LongitudinalGrid")
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

   // in
   MemBlock<double> memBlockIn(size, mem.get());
   // out and ref
   MemBlock<double> memBlockOut(size, mem.get());
   MemBlock<double> memBlockRef(size, mem.get());

   // views
   view_t in({memBlockIn.data(), memBlockIn.size()}, dimensions, pointers, indices);
   // out and ref
   view_t out({memBlockOut.data(), memBlockOut.size()}, dimensions, pointers, indices);
   view_t ref({memBlockRef.data(), memBlockRef.size()}, dimensions, pointers, indices);

   double scaling = 0.75;

   ::QuICC::Internal::Array igrid;
   ::QuICC::Internal::Array iweights;
   ::QuICC::Polynomial::Quadrature::LegendreRule quad;
   quad.computeQuadrature(igrid, iweights, out.dims()[1]);

   // init col 0, lay 0
   for (std::size_t m = 0; m < M; ++m)
   {
      auto mnl = m;
      in[mnl] = 1.0;
      ref[mnl] = scaling * (in[mnl] * std::sin(QuICC::Internal::cast(igrid[0])));
   }

   // init col 1, lay 2
   for (std::size_t m = 0; m < M; ++m)
   {
      auto mnl = m + M;
      in[mnl] = 1.0;
      ref[mnl] = scaling * (in[mnl] * std::sin(QuICC::Internal::cast(igrid[1])));
   }

   // const grid mul op
   using namespace QuICC::Slicewise::Cpu;
   using namespace QuICC::Slicewise;
   auto mulGridOp =
      std::make_unique<Op<1, ::QuICC::Polynomial::Quadrature::LegendreRule, MulSinFunctor<double>,
         view_t, view_t>>(MulSinFunctor<double>(scaling), mem);

   mulGridOp->apply(out, in);

   // check
   for (std::uint64_t i = 0; i < out.size(); ++i)
   {
      CHECK(out[i] == ref[i]);
   }
}
