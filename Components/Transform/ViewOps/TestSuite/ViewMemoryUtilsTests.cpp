#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "ViewOps/ViewMemoryUtils.hpp"

using namespace QuICC::Memory;
using namespace QuICC::View;
TEST_CASE("Do nothing", "[DoNothing]")
{
    // host mem block
    QuICC::Memory::Cpu::NewDelete mem;
    QuICC::Memory::MemBlock<double> memBlock(100, &mem);

    // view
    ViewBase<double> view(memBlock.data(), memBlock.size());

    // copy for reference
    ViewBase<double> originalView(view);

    // convert
    tempOnHostMemorySpace converter(view, TransferMode::read | TransferMode::block);

    // check
    CHECK(view.data() == originalView.data());
    CHECK(view.size() == originalView.size());
}


#ifdef QUICC_HAS_CUDA_BACKEND
TEST_CASE("Move read only", "[MoveRO]")
{
    // host mem block
    QuICC::Memory::Cuda::Malloc mem;
    QuICC::Memory::MemBlock<double> memBlock(100, &mem);

    // view
    ViewBase<double> view(memBlock.data(), memBlock.size());

    CHECK(QuICC::Cuda::isDeviceMemory(view.data()));

    // copy for reference
    ViewBase<double> originalView(view);

    {
        // convert
        tempOnHostMemorySpace converter(view, TransferMode::read);

        // check converted
        CHECK(view.data() != originalView.data());
        CHECK(!QuICC::Cuda::isDeviceMemory(view.data()));
        CHECK(view.size() == originalView.size());
    }

    // check restored
    CHECK(view.data() == originalView.data());
    CHECK(view.size() == originalView.size());
}

TEST_CASE("Move write and read", "[MoveWR]")
{
    // host mem block
    QuICC::Memory::Cuda::Malloc mem;
    QuICC::Memory::MemBlock<double> memBlock(100, &mem);

    // view
    ViewBase<double> view(memBlock.data(), memBlock.size());

    CHECK(QuICC::Cuda::isDeviceMemory(view.data()));

    // copy for reference
    ViewBase<double> originalView(view);

    {
        // convert
        tempOnHostMemorySpace converter(view, TransferMode::write | TransferMode::block);

        // check converted
        CHECK(view.data() != originalView.data());
        CHECK(!QuICC::Cuda::isDeviceMemory(view.data()));
        CHECK(view.size() == originalView.size());

        for (std::size_t i = 0; i < view.size(); ++i)
        {
            view[i] = i;
        }
    }

    // check restored
    CHECK(view.data() == originalView.data());
    CHECK(view.size() == originalView.size());

    {
        // convert
        tempOnHostMemorySpace converter(view, TransferMode::read | TransferMode::block);

        // check converted
        CHECK(view.data() != originalView.data());
        CHECK(!QuICC::Cuda::isDeviceMemory(view.data()));
        CHECK(view.size() == originalView.size());

        for (std::size_t i = 0; i < view.size(); ++i)
        {
            CHECK(view[i] == i);
        }
    }

    // check restored
    CHECK(view.data() == originalView.data());
    CHECK(view.size() == originalView.size());
}

#endif
