/**
 * @file Malloc.hpp
 * @brief
 */

#pragma once

// External includes
//
#include <utility>

// Project includes
//
#include "Memory/MemoryResource.hpp"

namespace QuICC {
namespace Memory {
/// @brief This namespace provides cuda only memory_resources
namespace Cuda {

/**
 *  @brief This class provides a concrete implementation of a memory_resources.
 * It is a wrapper to cuda operators cudaMalloc and cudaFree.
 *
 * */
class Malloc : public memory_resource
{
public:
    /// @brief ctors
    explicit Malloc() = default;
    /// @brief dtor
    ~Malloc() = default;
    /// @brief delete copy ctor
    /// @param
    Malloc (const Malloc&) = delete;
    /// @brief delete assignment operator
    /// @param
    /// @return
    Malloc& operator= (const Malloc&) = delete;
protected:
    /// @brief allocates memory
    /// @param bytes  size in bytes of the memory to be allocated
    /// @param alignment optional alignment in bytes
    /// @return pointer to allocated memory
    void* do_allocate(std::size_t bytes, std::size_t alignment) override;
    /// @brief deallocates memory
    /// @param p pointer to memory
    /// @param bytes  size in bytes of the memory to be allocated
    /// @param alignment optional alignment in bytes
    void do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) override;
    /// @brief compare for equality with another memory_resource
    /// @param other ref to memory resource
    /// @return true if equivalent
    bool do_is_equal(const QuICC::Memory::memory_resource& other) const noexcept override;
};

} // namespace Cuda
} // namespace Memory
} // namespace QuICC
