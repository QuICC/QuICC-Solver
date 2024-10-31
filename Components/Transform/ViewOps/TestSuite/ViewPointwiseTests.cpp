#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "ViewOps/Pointwise/Cpu/Pointwise.hpp"
#include "ViewOps/Pointwise/Functors.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"

using namespace QuICC::Memory;
using namespace QuICC::View;

TEST_CASE("Square Cpu Lambda", "SquareCpuLambda")
{
   // host mem block
   QuICC::Memory::Cpu::NewDelete mem;
   QuICC::Memory::MemBlock<double> memBlock(100, &mem);

   // view
   using view_t = ViewBase<double>;
   view_t view(memBlock.data(), memBlock.size());

   // init
   double val = 2.0;
   for (std::uint64_t i = 0; i < view.size(); ++i)
   {
      view[i] = val;
   }

   // square

   // cannot pass directly a lambda since we don't know the type
   auto squareLa = [](double in) { return in * in; };
   using namespace QuICC::Pointwise::Cpu;
   auto squareOp =
      std::make_unique<Op<decltype(squareLa), view_t, view_t>>(squareLa);

   squareOp->apply(view, view);

   // check
   for (std::uint64_t i = 0; i < view.size(); ++i)
   {
      CHECK(view[i] == val * val);
   }
}

TEST_CASE("Square Cpu Functor", "SquareCpuFunctor")
{
   // host mem block
   QuICC::Memory::Cpu::NewDelete mem;
   QuICC::Memory::MemBlock<double> memBlock(100, &mem);

   // view
   using view_t = ViewBase<double>;
   view_t view(memBlock.data(), memBlock.size());

   // init
   double val = 2.0;
   for (std::uint64_t i = 0; i < view.size(); ++i)
   {
      view[i] = val;
   }

   // square
   using namespace QuICC::Pointwise::Cpu;
   using namespace QuICC::Pointwise;
   auto squareOp = std::make_unique<Op<SquareFunctor<double>, view_t, view_t>>(
      SquareFunctor<double>());

   squareOp->apply(view, view);

   // check
   for (std::uint64_t i = 0; i < view.size(); ++i)
   {
      CHECK(view[i] == val * val);
   }
}

TEST_CASE("Abs2 Cpu Functor", "Abs2CpuFunctor")
{
   auto constexpr S = 5;
   // host mem block
   QuICC::Memory::Cpu::NewDelete mem;
   QuICC::Memory::MemBlock<std::complex<double>> memBlockIn(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockOut(S, &mem);

   // view in
   using viewIn_t = ViewBase<std::complex<double>>;
   viewIn_t viewIn(memBlockIn.data(), memBlockIn.size());

   // view out
   using viewOut_t = ViewBase<double>;
   viewOut_t viewOut(memBlockOut.data(), memBlockOut.size());

   // init
   std::complex<double> val{2.0, -2.0};
   for (std::uint64_t i = 0; i < viewIn.size(); ++i)
   {
      viewIn[i] = val;
   }

   // abs2
   using namespace QuICC::Pointwise::Cpu;
   using namespace QuICC::Pointwise;
   auto abs2Op = std::make_unique<Op<Abs2Functor<double>, viewOut_t, viewIn_t>>(
      Abs2Functor<double>());

   abs2Op->apply(viewOut, viewIn);

   // check
   for (std::uint64_t i = 0; i < viewOut.size(); ++i)
   {
      CHECK(viewOut[i] == val.real() * val.real() + val.imag() * val.imag());
   }
}

TEST_CASE("Add Cpu Functor", "AddCpuFunctor")
{
   auto constexpr S = 5;
   // host mem block
   QuICC::Memory::Cpu::NewDelete mem;
   QuICC::Memory::MemBlock<double> memBlockIn0(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockIn1(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockOut(S, &mem);

   // view in
   using viewIn_t = ViewBase<double>;
   viewIn_t viewIn0(memBlockIn0.data(), memBlockIn0.size());
   viewIn_t viewIn1(memBlockIn1.data(), memBlockIn1.size());

   // view out
   using viewOut_t = ViewBase<double>;
   viewOut_t viewOut(memBlockOut.data(), memBlockOut.size());

   // init
   double val0 = 2.0;
   double val1 = 3.0;
   for (std::uint64_t i = 0; i < S; ++i)
   {
      viewIn0[i] = val0;
      viewIn1[i] = val1;
   }

   // add
   using namespace QuICC::Pointwise::Cpu;
   using namespace QuICC::Pointwise;
   auto addOp =
      std::make_unique<Op<AddFunctor<double>, viewOut_t, viewIn_t, viewIn_t>>(
         AddFunctor<double>());

   addOp->apply(viewOut, viewIn0, viewIn1);

   // check
   for (std::uint64_t i = 0; i < S; ++i)
   {
      CHECK(viewOut[i] == val0 + val1);
   }
}

TEST_CASE("Cross Cpu Functor", "CrossCpuFunctor")
{
   auto constexpr S = 5;
   // host mem block
   QuICC::Memory::Cpu::NewDelete mem;
   QuICC::Memory::MemBlock<double> memBlockU0(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockU1(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockU2(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockV0(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockV1(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockV2(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockC0(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockC1(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockC2(S, &mem);

   // view in
   using view_t = ViewBase<double>;
   view_t viewU0(memBlockU0.data(), memBlockU0.size());
   view_t viewU1(memBlockU1.data(), memBlockU1.size());
   view_t viewU2(memBlockU2.data(), memBlockU2.size());
   view_t viewV0(memBlockV0.data(), memBlockV0.size());
   view_t viewV1(memBlockV1.data(), memBlockV1.size());
   view_t viewV2(memBlockV2.data(), memBlockV2.size());

   // view Cross
   view_t viewC0(memBlockC0.data(), memBlockC0.size());
   view_t viewC1(memBlockC1.data(), memBlockC1.size());
   view_t viewC2(memBlockC2.data(), memBlockC2.size());

   // init
   for (std::uint64_t i = 0; i < S; ++i)
   {
      viewU0[i] = 1.6;
      viewU1[i] = 2.5;
      viewU2[i] = 3.4;
      viewV0[i] = 4.3;
      viewV1[i] = 5.2;
      viewV2[i] = 6.1;
   }

   // Cross
   double scaling = 0.25;
   using namespace QuICC::Pointwise::Cpu;
   using namespace QuICC::Pointwise;
   auto addOp =
      std::make_unique<Op<CrossCompFunctor<double>, view_t, view_t, view_t, view_t, view_t>>(
         CrossCompFunctor<double>(scaling));

   addOp->apply(viewC0, viewU1, viewU2, viewV1, viewV2);

   // check
   auto eps = 10*std::numeric_limits<double>::epsilon();
   for (std::uint64_t i = 0; i < S; ++i)
   {
      CHECK(std::abs(viewC0[i] - scaling * (viewU1[i] * viewV2[i] - viewU2[i] * viewV1[i])) < eps);
   }
}

TEST_CASE("Dot Cpu Functor", "DotCpuFunctor")
{
   auto constexpr S = 5;
   // host mem block
   QuICC::Memory::Cpu::NewDelete mem;
   QuICC::Memory::MemBlock<double> memBlockU0(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockU1(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockU2(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockV0(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockV1(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockV2(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockOut(S, &mem);

   // view in
   using view_t = ViewBase<double>;
   view_t viewU0(memBlockU0.data(), memBlockU0.size());
   view_t viewU1(memBlockU1.data(), memBlockU1.size());
   view_t viewU2(memBlockU2.data(), memBlockU2.size());
   view_t viewV0(memBlockV0.data(), memBlockV0.size());
   view_t viewV1(memBlockV1.data(), memBlockV1.size());
   view_t viewV2(memBlockV2.data(), memBlockV2.size());

   // view Dot
   view_t viewOut(memBlockOut.data(), memBlockOut.size());
   // init
   for (std::uint64_t i = 0; i < S; ++i)
   {
      viewU0[i] = 1.6;
      viewU1[i] = 2.5;
      viewU2[i] = 3.4;
      viewV0[i] = 4.3;
      viewV1[i] = 5.2;
      viewV2[i] = 6.1;
   }

   // Dot
   double scaling = 0.25;
   using namespace QuICC::Pointwise::Cpu;
   using namespace QuICC::Pointwise;
   auto addOp =
      std::make_unique<Op<DotFunctor<double>, view_t, view_t, view_t, view_t, view_t, view_t, view_t>>(
         DotFunctor<double>(scaling));

   addOp->apply(viewOut, viewU0, viewU1, viewU2, viewV0, viewV1, viewV2);

   // check
   auto eps = 10*std::numeric_limits<double>::epsilon();
   for (std::uint64_t i = 0; i < S; ++i)
   {
      CHECK(std::abs(viewOut[i] - scaling * (viewU0[i] * viewV0[i]
      + viewU1[i] * viewV1[i] + viewU2[i] * viewV2[i])) < eps);
   }
}