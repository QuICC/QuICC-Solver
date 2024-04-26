#include <catch2/catch.hpp>

#include "ViewOps/Pointwise/Cpu/Pointwise.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#include "ViewOps/Pointwise/Cuda/Pointwise.hpp"
#include "ViewOps/Pointwise/Functors.hpp"

using namespace QuICC::Memory;
using namespace QuICC::View;


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
   auto squareOp = std::make_unique<Op<SquareFunctor<double>, view_t, view_t>>(
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
    using inTy = std::complex<double>;
    auto constexpr S = 5;
    // host mem block
    QuICC::Memory::Cpu::NewDelete mem;
    QuICC::Memory::MemBlock<inTy> memBlockIn(S, &mem);
    QuICC::Memory::MemBlock<double> memBlockOut(S, &mem);

    // device mem block
    QuICC::Memory::Cuda::Malloc memDev;
    QuICC::Memory::MemBlock<inTy> memBlockInDev(S, &memDev);
    QuICC::Memory::MemBlock<double> memBlockOutDev(S, &memDev);

    // view in
    using viewIn_t = ViewBase<inTy>;
    viewIn_t viewIn(memBlockIn.data(), memBlockIn.size());
    viewIn_t viewInDev(memBlockInDev.data(), memBlockInDev.size());

    // view out
    using viewOut_t = ViewBase<double>;
    viewOut_t viewOut(memBlockOut.data(), memBlockOut.size());
    viewOut_t viewOutDev(memBlockOutDev.data(), memBlockOutDev.size());

    // init
    inTy val{2.0, -2.0};
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
    auto abs2Op = std::make_unique<Op<Abs2Functor<double>, viewOut_t, viewIn_t>>(
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

TEST_CASE("Add complex Gpu", "AddComplexGpu")
{
    // c = a + b
    auto constexpr S = 5;
    // host mem block
    QuICC::Memory::Cpu::NewDelete mem;
    using Ty = cuda::std::complex<double>;
    QuICC::Memory::MemBlock<Ty> memBlockA(S, &mem);
    QuICC::Memory::MemBlock<Ty> memBlockB(S, &mem);
    QuICC::Memory::MemBlock<Ty> memBlockC(S, &mem);

    // device mem block
    QuICC::Memory::Cuda::Malloc memDev;
    QuICC::Memory::MemBlock<Ty> memBlockADev(S, &memDev);
    QuICC::Memory::MemBlock<Ty> memBlockBDev(S, &memDev);
    QuICC::Memory::MemBlock<Ty> memBlockCDev(S, &memDev);

    // host view
    using view_t = ViewBase<Ty>;
    view_t viewA(memBlockA.data(), memBlockA.size());
    view_t viewB(memBlockB.data(), memBlockB.size());
    view_t viewC(memBlockC.data(), memBlockC.size());

    // device view
    view_t viewADev(memBlockADev.data(), memBlockADev.size());
    view_t viewBDev(memBlockBDev.data(), memBlockBDev.size());
    view_t viewCDev(memBlockCDev.data(), memBlockCDev.size());

    // init
    Ty valA = {2.0, -1.0};
    Ty valB = {3.0, 5.0};
    for (std::uint64_t i = 0; i < viewA.size(); ++i)
    {
        viewA[i] = valA;
        viewB[i] = valB;
    }

    // Copy to gpu
    cudaErrChk(cudaMemcpy(viewADev.data(), viewA.data(),
        S * sizeof(typename view_t::ScalarType), cudaMemcpyHostToDevice));
    cudaErrChk(cudaMemcpy(viewBDev.data(), viewB.data(),
        S * sizeof(typename view_t::ScalarType), cudaMemcpyHostToDevice));

    // Add

    // we need to define a functor instead of using lambda to explicitly
    // instantiate backend
    using namespace QuICC::Pointwise::Cuda;
    using namespace QuICC::Pointwise;
    auto addOp = std::make_unique<Op<AddFunctor<Ty>, view_t, view_t, view_t>>(
        AddFunctor<Ty>());

    addOp->apply(viewCDev, viewADev, viewBDev);

    // copy back
    cudaErrChk(cudaMemcpy(viewC.data(), viewCDev.data(),
        S * sizeof(typename view_t::ScalarType), cudaMemcpyDeviceToHost));

    // check
    for (std::uint64_t i = 0; i < viewC.size(); ++i)
    {
        CHECK(viewC[i].real() == (valA + valB).real());
        CHECK(viewC[i].imag() == (valA + valB).imag());
    }
}

TEST_CASE("Sub complex Gpu", "SubComplexGpu")
{
    // c = a - b
    auto constexpr S = 5;
    // host mem block
    QuICC::Memory::Cpu::NewDelete mem;
    using Ty = cuda::std::complex<double>;
    QuICC::Memory::MemBlock<Ty> memBlockA(S, &mem);
    QuICC::Memory::MemBlock<Ty> memBlockB(S, &mem);
    QuICC::Memory::MemBlock<Ty> memBlockC(S, &mem);

    // device mem block
    QuICC::Memory::Cuda::Malloc memDev;
    QuICC::Memory::MemBlock<Ty> memBlockADev(S, &memDev);
    QuICC::Memory::MemBlock<Ty> memBlockBDev(S, &memDev);
    QuICC::Memory::MemBlock<Ty> memBlockCDev(S, &memDev);

    // host view
    using view_t = ViewBase<Ty>;
    view_t viewA(memBlockA.data(), memBlockA.size());
    view_t viewB(memBlockB.data(), memBlockB.size());
    view_t viewC(memBlockC.data(), memBlockC.size());

    // device view
    view_t viewADev(memBlockADev.data(), memBlockADev.size());
    view_t viewBDev(memBlockBDev.data(), memBlockBDev.size());
    view_t viewCDev(memBlockCDev.data(), memBlockCDev.size());

    // init
    Ty valA = {2.0, -1.0};
    Ty valB = {3.0, 5.0};
    for (std::uint64_t i = 0; i < viewA.size(); ++i)
    {
        viewA[i] = valA;
        viewB[i] = valB;
    }

    // Copy to gpu
    cudaErrChk(cudaMemcpy(viewADev.data(), viewA.data(),
        S * sizeof(typename view_t::ScalarType), cudaMemcpyHostToDevice));
    cudaErrChk(cudaMemcpy(viewBDev.data(), viewB.data(),
        S * sizeof(typename view_t::ScalarType), cudaMemcpyHostToDevice));

    // Sub

    // we need to define a functor instead of using lambda to explicitly
    // instantiate backend
    using namespace QuICC::Pointwise::Cuda;
    using namespace QuICC::Pointwise;
    auto subOp = std::make_unique<Op<SubFunctor<Ty>, view_t, view_t, view_t>>(
        SubFunctor<Ty>());

    subOp->apply(viewCDev, viewADev, viewBDev);

    // copy back
    cudaErrChk(cudaMemcpy(viewC.data(), viewCDev.data(),
        S * sizeof(typename view_t::ScalarType), cudaMemcpyDeviceToHost));

    // check
    for (std::uint64_t i = 0; i < viewC.size(); ++i)
    {
        CHECK(viewC[i].real() == (valA - valB).real());
        CHECK(viewC[i].imag() == (valA - valB).imag());
    }
}

