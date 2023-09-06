#include <catch2/catch.hpp>

#include "Memory/Memory.hpp"
#include "Memory/Cuda/Malloc.hpp"
#include "Cuda/CudaUtil.hpp"

__global__
void set1(double* ptr, std::size_t size)
{
    const std::size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index < size) ptr[index] = 1.0;
}


TEST_CASE("MemoryBlock using CudaMalloc", "[MemoryBlockCudaMalloc]")
{
    constexpr std::size_t S = 10;
    std::array<double, S> block_h;
    QuICC::Memory::Cuda::Malloc mem_res;
    QuICC::Memory::MemBlock<std::byte> block_d(S*sizeof(double), &mem_res);
    
    // set to ones
    const unsigned int blockSize = 32;
    const unsigned int numBlocks = (S + blockSize - 1) / blockSize;
    set1<<<blockSize, numBlocks>>>(reinterpret_cast<double*>(block_d.data()), S);
    cudaErrChk(cudaDeviceSynchronize());

    // copy back 
    cudaErrChk(cudaMemcpy(block_h.data(), block_d.data(),
        S*sizeof(double), cudaMemcpyDeviceToHost));
    
    // check
    for (std::size_t i = 0; i < S; ++i)
    {
        CHECK(block_h[i] == 1.0);
    }

}

TEST_CASE("MemoryBlock using CudaMalloc and Move", "[MemoryBlockCudaMallocMove]")
{
    constexpr std::size_t S = 10;
    std::array<double, S> block_h;
    QuICC::Memory::Cuda::Malloc mem_res;
    QuICC::Memory::MemBlock<std::byte> block_d;
    QuICC::Memory::MemBlock<std::byte> tmp(S*sizeof(double), &mem_res);
    
    block_d = std::move(tmp);

    // set to ones
    const unsigned int blockSize = 32;
    const unsigned int numBlocks = (S + blockSize - 1) / blockSize;
    set1<<<blockSize, numBlocks>>>(reinterpret_cast<double*>(block_d.data()), S);
    cudaErrChk(cudaDeviceSynchronize());

    // copy back 
    cudaErrChk(cudaMemcpy(block_h.data(), block_d.data(),
        S*sizeof(double), cudaMemcpyDeviceToHost));
    
    // check
    for (std::size_t i = 0; i < S; ++i)
    {
        CHECK(block_h[i] == 1.0);
    }

}
