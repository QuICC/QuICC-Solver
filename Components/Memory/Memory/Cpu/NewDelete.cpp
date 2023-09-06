/**
 * @file NewDelete.cpp
 * @brief
 */


// External includes
//
#include <cstddef>
#ifndef NDEBUG
#include <iostream>
#endif
#include <new>

// Project includes
//
#include "NewDelete.hpp"

namespace QuICC {
namespace Memory {
namespace Cpu {

void* NewDelete::do_allocate(std::size_t bytes, std::size_t alignment)
{
    void* ptr{nullptr};

    // new
    ptr = ::operator new(bytes, static_cast<std::align_val_t>(alignment));

    #ifndef NDEBUG
    std::cout << "New, bytes: " << bytes << '\n';
    #endif

    return ptr;
}

void NewDelete::do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment)
{
    if(ptr != nullptr)
    {
        #ifndef NDEBUG
        std::cout << "Delete, bytes: " << bytes << '\n';
        #endif
        ::operator delete(ptr, bytes, static_cast<std::align_val_t>(alignment));
    }
}

bool NewDelete::do_is_equal(const QuICC::Memory::memory_resource& other) const noexcept
{
    if( this == &other ) return true;
    const auto* const op = dynamic_cast<const NewDelete*>(&other);
    return op != nullptr;
}

} // namespace Cpu
} // namespace Memory
} // namespace QuICC

