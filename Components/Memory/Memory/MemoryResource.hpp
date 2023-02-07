/**
 * @file MemoryResource.hpp
 * @brief mimics https://en.cppreference.com/w/cpp/memory/memory_resource 
 * designed to work for both device and host code
 */

#pragma once

// External includes
//
#include <cstddef>

// Project includes
//

namespace QuICC {
namespace Memory {

/** @brief This class is an abstract interface to an unbounded set of classes
 * encapsulating memory resources.
 * It is used instead of a std::memory_resource so that it can be extended to be used in device functions.
 * https://en.cppreference.com/w/cpp/memory/memory_resource
 */
class memory_resource
{
public:
    /// @brief destruct a memory_resource
    virtual ~memory_resource() = default;

    /// @brief allocates memory
    /// @param bytes  size in bytes of the memory to be allocated
    /// @param alignment optional alignment in bytes
    /// @return pointer to allocated memory
    void* allocate( std::size_t bytes,
        std::size_t alignment = alignof(std::max_align_t) )
    {
        return do_allocate(bytes, alignment);
    }

    /// @brief deallocates memory
    /// @param p pointer to memory
    /// @param bytes  size in bytes of the memory to be allocated
    /// @param alignment optional alignment in bytes
    void deallocate( void* p, std::size_t bytes,
        std::size_t alignment = alignof(std::max_align_t) )
    {
        do_deallocate(p, bytes, alignment);
    }

    /// @brief compare for equality with another memory_resource
    /// @param other ref to memory resource
    /// @return true if equivalent
    bool is_equal( const memory_resource& other ) const noexcept
    {
        return do_is_equal(other);
    }

private:

    /// @brief allocates memory, overwritten by derived class
    /// @param bytes  size in bytes of the memory to be allocated
    /// @param alignment optional alignment in bytes
    /// @return pointer to allocated memory
    virtual void* do_allocate( std::size_t bytes, std::size_t alignment ) = 0;

    /// @brief deallocates memory, overwritten by derived class
    /// @param p pointer to memory
    /// @param bytes  size in bytes of the memory to be allocated
    /// @param alignment optional alignment in bytes
    virtual void do_deallocate( void* p, std::size_t bytes, std::size_t alignment ) = 0;

    /// @brief compare for equality with another memory_resource, overwritten by derived class
    /// @param other ref to memory resource
    /// @return true if equivalent
    virtual bool do_is_equal( const memory_resource& other ) const noexcept = 0;
};

} // namespace Memory
} // namespace QuICC

