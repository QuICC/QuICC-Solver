#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "Memory/Memory.hpp"
#include "Memory/Cpu/NewDelete.hpp"
#include "Memory/Pensieve.hpp"

using namespace QuICC::Memory;

TEST_CASE("MemoryBlock", "[MemoryBlock]")
{
    constexpr std::size_t S = 10;
    Cpu::NewDelete mem_res;
    MemBlock<double> block_h(S, &mem_res);

    // set to ones
    for (std::size_t i = 0; i < S; ++i)
    {
        block_h.data()[i] = 1.0;
    }

    // check
    for (std::size_t i = 0; i < S; ++i)
    {
        CHECK(block_h.data()[i] == 1.0);
    }
}

TEST_CASE("MemoryBlockByte", "[MemoryBlockByte]")
{
    constexpr std::size_t S = 10;
    Cpu::NewDelete mem_res;
    MemBlock<std::byte> block_h(S*sizeof(double), &mem_res);

    // set to ones
    for (std::size_t i = 0; i < S; ++i)
    {
        reinterpret_cast<double*>(block_h.data())[i] = 1.0;
    }

    // check
    for (std::size_t i = 0; i < S; ++i)
    {
        CHECK(reinterpret_cast<double*>(block_h.data())[i] == 1.0);
    }
}

TEST_CASE("MemoryBlockMove", "[MemoryBlockMove]")
{
    constexpr std::size_t S = 10;
    Cpu::NewDelete mem_res;
    MemBlock<double> block_h(S, &mem_res);

    MemBlock<double> block_h2;

    CHECK(block_h2.data() == nullptr);
    CHECK(block_h2.size() == 0);

    block_h2 = std::move(block_h);

    CHECK(block_h.data() == nullptr);
    CHECK(block_h.size() == 0);

    CHECK(block_h2.data() != nullptr);
    CHECK(block_h2.size() == S);
}

TEST_CASE("Pensieve", "[Pensieve]")
{
    [[maybe_unused]] auto& mem = Pensieve<Cpu::NewDelete>::getInstance().getMem();
    // this should trigger a static assert
    // [[maybe_unused]] auto& mem2 = Pensieve<double>::getInstance().getMem();
}
