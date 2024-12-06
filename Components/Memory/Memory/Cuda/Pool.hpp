/**
 * @file Pool.hpp
 * @brief
 */

#pragma once

// External includes
//
#include <utility>
#include <vector>
#include <stack>
#include <cstdint>

// Project includes
//
#include "Memory/MemoryResource.hpp"

namespace QuICC {
namespace Memory {
/// @brief This namespace provides Cuda only memory_resources
namespace Cuda {


namespace details
{
    /// @brief block of the memory pool
    struct block
    {
        /// @brief pointer to memory block
        /// 8 bytes
        void* ptr = nullptr;
        /// @brief block size in bytes
        /// 4 bytes
        std::uint32_t blockSize = 0;
        /// @brief flag to check if in use
        /// 1 byte
        bool isValid = true;
    };

} // namespace details


/**
 *  @brief This class provides a concrete implementation of a memory_resources.
 * It is a simple memory pool implementation.
 * Alignement is fixed to cache line size.
 * The block size is grown dynamically to fit largest allocation.
 *
 * */
class Pool : public memory_resource
{
public:
    /// @brief ctor
    explicit Pool() = default;
    /// @brief dtor
    ~Pool();
    /// @brief delete copy ctor
    /// @param
    Pool (const Pool&) = delete;
    /// @brief delete assignment operator
    /// @param
    /// @return
    Pool& operator= (const Pool&) = delete;
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

    /// @brief memory blocks storage
    std::vector<details::block> _blocks;
    /// @brief stack keeping track of free blocks indeces
    std::stack<std::uint32_t> _freeBlocks;
    /// @brief alignement in bytes, cache line
    const std::uint32_t _alignment = 64;
    /// @brief current biggest block size
    std::size_t _maxBlockSize = 128;

};

} // namespace Cuda
} // namespace Memory
} // namespace QuICC

