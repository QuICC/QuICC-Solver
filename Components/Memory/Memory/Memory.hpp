/**
 * @file Memory.hpp
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
/// @brief This namespace provides methods and classes for memory management
namespace Memory {

/** @brief This class represents a single, contigous block of memory.
 *  MemBlock is the owner of the memory, it is not device aware and through a pmr style approach can implement a memory pool.
 *  It is used instead of a std container so that it can be extended to be used in device functions.
 *  @tparam T
 */
template<class T>
class MemBlock
{
public:
    /// @brief typedef to retrieve element data type
    using type = T;
    /// @brief default ctor (empty block)
    MemBlock() = default;
    /// @brief delete copy ctor
    MemBlock(const MemBlock&) = delete;
    /// @brief move ctor
    MemBlock(MemBlock&& old) :
        _data(std::exchange(old._data, nullptr)),
        _size(std::exchange(old._size, 0)),
        _mem_res(std::exchange(old._mem_res, nullptr)) {}

    /// @brief ctor
    /// @param size number of elements of type type
    /// @param mem memory resource pointer
    explicit MemBlock(std::size_t size, memory_resource* mem) :
        _mem_res(mem)
    {
        _data = reinterpret_cast<T*>(_mem_res->allocate(size*sizeof(T)));
        _size = size;
    }

    /// @brief dtor, release memory resource
    ~MemBlock()
    {
        _release();
    }
    
    /// @brief delete copy assignment
    MemBlock& operator=(const MemBlock&) = delete;
    /// @brief move assignement
    MemBlock& operator=(MemBlock&& old)
    {
        _release();
        _data = std::exchange(old._data, nullptr);
        _size = std::exchange(old._size, 0);
        _mem_res = std::exchange(old._mem_res, nullptr);
        return *this;
    }

    /// @brief access raw data pointer with const qualifier
    /// @return _data
    const T* operator*() const {return _data;}

    /// @brief access raw data pointer
    /// @return _data
    T* operator*() {return _data;}

    /// @brief access raw data pointer with const qualifier
    /// @return _data
    const T* data() const {return _data;}

    /// @brief access raw data pointer
    /// @return _data
    T* data() {return _data;}

    /// @brief get size in number of elements
    /// @return _size
    std::size_t size() const {return _size;}

private:

    /// @brief release resource if needed
    void _release()
    {
        if (_data != nullptr)
        {
            _mem_res->deallocate(_data, _size*sizeof(T));
            _size = 0;
            _data = nullptr;
        }
    }

    /// @brief pointer to raw data
    T* _data{nullptr};

    /// @brief size of block in number of elements
    std::size_t _size{0};

    /// @brief pointer to memory resource
    memory_resource* _mem_res{nullptr};
};

} // namespace Memory
} // namespace QuICC

