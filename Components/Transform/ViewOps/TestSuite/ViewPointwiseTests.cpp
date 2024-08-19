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
