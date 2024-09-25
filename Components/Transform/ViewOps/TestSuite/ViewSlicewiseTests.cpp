#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "ViewOps/Slicewise/Op.hpp"
#include "ViewOps/Slicewise/Functors.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"

using namespace QuICC::Memory;
using namespace QuICC::View;

TEST_CASE("Spherical Heat Advection", "SphericalHeatAdvection")
{
   // host mem block
   QuICC::Memory::Cpu::NewDelete mem;
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

   // vel
   MemBlock<double> memBlockVR(size, &mem);
   MemBlock<double> memBlockVTheta(size, &mem);
   MemBlock<double> memBlockVPhi(size, &mem);
   // gradT
   MemBlock<double> memBlockTdR(size, &mem);
   MemBlock<double> memBlockTdTheta(size, &mem);
   MemBlock<double> memBlockTdPhi(size, &mem);
   // out and ref
   MemBlock<double> memBlockOut(size, &mem);
   MemBlock<double> memBlockRef(size, &mem);
   // grid
   MemBlock<double> memBlockGrid(L, &mem);

   // views
   // vel
   view_t vR({memBlockVR.data(), memBlockVR.size()}, dimensions, pointers, indices);
   view_t vTheta({memBlockVTheta.data(), memBlockVTheta.size()}, dimensions, pointers, indices);
   view_t vPhi({memBlockVPhi.data(), memBlockVPhi.size()}, dimensions, pointers, indices);
   // gradT
   view_t TdR({memBlockTdR.data(), memBlockTdR.size()}, dimensions, pointers, indices);
   view_t TdTheta({memBlockTdTheta.data(), memBlockTdTheta.size()}, dimensions, pointers, indices);
   view_t TdPhi({memBlockTdPhi.data(), memBlockTdPhi.size()}, dimensions, pointers, indices);
   // out and ref
   view_t out({memBlockOut.data(), memBlockOut.size()}, dimensions, pointers, indices);
   view_t ref({memBlockRef.data(), memBlockRef.size()}, dimensions, pointers, indices);
   // grid
   ViewBase<double> grid(memBlockGrid.data(), memBlockGrid.size());

   // init grid
   grid[0] = 0.1;
   grid[1] = 0.5;
   grid[2] = 0.9;

   double scaling = 0.7;

   // init col 0, lay 0
   for (std::size_t m = 0; m < M; ++m)
   {
      auto mnl = m;
      vR[mnl] = 1.0;
      vTheta[mnl] = 2.0;
      vPhi[mnl] = 3.0;
      TdR[mnl] = 4.0;
      TdTheta[mnl] = -5.0;
      TdPhi[mnl] = 6.0;
      ref[mnl] = scaling * (vR[mnl]*TdR[mnl] + vTheta[mnl]*TdTheta[mnl] + vPhi[mnl]*TdPhi[mnl]
         - grid[0] * vR[mnl]);
   }

   // init col 1, lay 2
   for (std::size_t m = 0; m < M; ++m)
   {
      auto mnl = m + M;
      vR[mnl] = 1.0;
      vTheta[mnl] = 2.0;
      vPhi[mnl] = 3.0;
      TdR[mnl] = 4.0;
      TdTheta[mnl] = 5.0;
      TdPhi[mnl] = 6.0;
      ref[mnl] = scaling * (vR[mnl]*TdR[mnl] + vTheta[mnl]*TdTheta[mnl] + vPhi[mnl]*TdPhi[mnl]
         - grid[2] * vR[mnl]);
   }

   // advection op
   using namespace QuICC::Slicewise::Cpu;
   using namespace QuICC::Slicewise;
   auto advOp =
      std::make_unique<Op<SphericalHeatAdvection<double>,
         view_t, view_t, view_t, view_t, view_t, view_t, view_t>>(
         SphericalHeatAdvection<double>(grid, scaling));

   advOp->apply(out, vR, vTheta, vPhi, TdR, TdTheta, TdPhi);

   // check
   for (std::uint64_t i = 0; i < out.size(); ++i)
   {
      CHECK(out[i] == ref[i]);
   }
}
