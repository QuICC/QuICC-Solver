/**
 * @file Malloc.cu
 * @brief
 */


// External includes
//
#ifndef NDEBUG
#include <iostream>
#endif

// Project includes
//
#include "Malloc.hpp"
#include "Cuda/CudaUtil.hpp"

namespace QuICC {
namespace Memory {
namespace Cuda {

void* Malloc::do_allocate(std::size_t bytes, std::size_t alignment)
{
    void* ptr{nullptr};
    cudaErrChk(cudaMalloc(reinterpret_cast<void**>(&ptr), bytes));
    #ifndef NDEBUG
    std::cout << "cudaMalloc, bytes: " << bytes << '\n';
    #endif

    return ptr;
}

void Malloc::do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment)
{
    if(ptr != nullptr)
    {    
        #ifndef NDEBUG
        std::cout << "cudaFree, bytes: " << bytes << '\n';
        #endif
        cudaErrChk(cudaFree(ptr));
        ptr = nullptr;
    }
}

bool Malloc::do_is_equal(const QuICC::Memory::memory_resource& other) const noexcept
{
    if( this == &other ) return true;
    const auto* const op = dynamic_cast<const Malloc*>(&other);
    return op != nullptr;
}

} // namespace Cuda
} // namespace Memory
} // namespace QuICC

