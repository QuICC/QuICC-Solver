#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "Memory/Memory.hpp"
#include "Memory/Cpu/NewDelete.hpp"

TEST_CASE("MemoryBlock", "[MemoryBlock]")
{
    constexpr std::size_t S = 10;
    QuICC::Memory::Cpu::NewDelete mem_res;
    QuICC::Memory::MemBlock<double> block_h(S, &mem_res);
    
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
    QuICC::Memory::Cpu::NewDelete mem_res;
    QuICC::Memory::MemBlock<std::byte> block_h(S*sizeof(double), &mem_res);
    
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
