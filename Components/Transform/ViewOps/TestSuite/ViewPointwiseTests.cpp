#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "ViewOps/Pointwise/Cpu/Pointwise.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "ViewOps/Pointwise/Cuda/Pointwise.hpp"
#endif
#include "ViewOps/Pointwise/Functors.hpp"

using namespace QuICC::Memory;
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
      std::make_unique<Op<view_t, view_t, decltype(squareLa)>>(squareLa);

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
   auto squareOp = std::make_unique<Op<view_t, view_t, SquareFunctor<double>>>(
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
   auto abs2Op = std::make_unique<Op<viewOut_t, viewIn_t, Abs2Functor<double>>>(
      Abs2Functor<double>());

   abs2Op->apply(viewOut, viewIn);

   // check
   for (std::uint64_t i = 0; i < viewOut.size(); ++i)
   {
      CHECK(viewOut[i] == val.real() * val.real() + val.imag() * val.imag());
   }
}

#ifdef QUICC_HAS_CUDA_BACKEND

TEST_CASE("Square Gpu", "SquareGpu")
{
   auto constexpr S = 5;
   // host mem block
   QuICC::Memory::Cpu::NewDelete mem;
   QuICC::Memory::MemBlock<double> memBlock(S, &mem);

   // device mem block
   QuICC::Memory::Cuda::Malloc memDev;
   QuICC::Memory::MemBlock<double> memBlockDev(S, &memDev);

   // host view
   using view_t = ViewBase<double>;
   view_t view(memBlock.data(), memBlock.size());

   // device view
   view_t viewDev(memBlockDev.data(), memBlockDev.size());

   // init
   double val = 2.0;
   for (std::uint64_t i = 0; i < view.size(); ++i)
   {
      view[i] = val;
   }

   // Copy to gpu
   cudaErrChk(cudaMemcpy(viewDev.data(), view.data(),
      S * sizeof(typename view_t::ScalarType), cudaMemcpyHostToDevice));

   // square

   // we need to define a functor instead of using lambda to explicitly
   // instantiate backend
   using namespace QuICC::Pointwise::Cuda;
   using namespace QuICC::Pointwise;
   auto squareOp = std::make_unique<Op<view_t, view_t, SquareFunctor<double>>>(
      SquareFunctor<double>());

   squareOp->apply(viewDev, viewDev);

   // copy back
   cudaErrChk(cudaMemcpy(view.data(), viewDev.data(),
      S * sizeof(typename view_t::ScalarType), cudaMemcpyDeviceToHost));

   // check
   for (std::uint64_t i = 0; i < view.size(); ++i)
   {
      CHECK(view[i] == val * val);
   }
}

TEST_CASE("Abs2 Gpu", "Abs2Gpu")
{
   auto constexpr S = 5;
   // host mem block
   QuICC::Memory::Cpu::NewDelete mem;
   QuICC::Memory::MemBlock<std::complex<double>> memBlockIn(S, &mem);
   QuICC::Memory::MemBlock<double> memBlockOut(S, &mem);

   // device mem block
   QuICC::Memory::Cuda::Malloc memDev;
   QuICC::Memory::MemBlock<std::complex<double>> memBlockInDev(S, &memDev);
   QuICC::Memory::MemBlock<double> memBlockOutDev(S, &memDev);

   // view in
   using viewIn_t = ViewBase<std::complex<double>>;
   viewIn_t viewIn(memBlockIn.data(), memBlockIn.size());
   viewIn_t viewInDev(memBlockInDev.data(), memBlockInDev.size());

   // view out
   using viewOut_t = ViewBase<double>;
   viewOut_t viewOut(memBlockOut.data(), memBlockOut.size());
   viewOut_t viewOutDev(memBlockOutDev.data(), memBlockOutDev.size());

   // init
   std::complex<double> val{2.0, -2.0};
   for (std::uint64_t i = 0; i < viewIn.size(); ++i)
   {
      viewIn[i] = val;
   }

   // Copy to gpu
   cudaErrChk(cudaMemcpy(viewInDev.data(), viewIn.data(),
      S * sizeof(typename viewIn_t::ScalarType), cudaMemcpyHostToDevice));

   // abs2

   // we need to define a functor instead of using lambda to explicitly
   // instantiate backend
   using namespace QuICC::Pointwise::Cuda;
   using namespace QuICC::Pointwise;
   auto abs2Op = std::make_unique<Op<viewOut_t, viewIn_t, Abs2Functor<double>>>(
      Abs2Functor<double>());

   abs2Op->apply(viewOutDev, viewInDev);

   // copy back
   cudaErrChk(cudaMemcpy(viewOut.data(), viewOutDev.data(),
      S * sizeof(typename viewOut_t::ScalarType), cudaMemcpyDeviceToHost));

   // check
   for (std::uint64_t i = 0; i < viewOut.size(); ++i)
   {
      CHECK(viewOut[i] == val.real() * val.real() + val.imag() * val.imag());
   }
}


#endif
