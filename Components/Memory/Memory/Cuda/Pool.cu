/**
 * @file Pool.cpp
 * @brief
 */

// External includes
//
#include <cstddef>
#ifndef NDEBUG
#include <iostream>
#endif
#include <cassert>

// Project includes
//
#include "Memory/Cuda/Pool.hpp"
#include "Cuda/CudaUtil.hpp"

namespace QuICC {
namespace Memory {
namespace Cuda {

Pool::~Pool()
{
    // deallocate all blocks
    for(std::uint32_t i = 0; i < _blocks.size(); ++i)
    {
        if(_blocks[i].ptr != nullptr)
        {
            // ::operator delete(_blocks[i].ptr, _blocks[i].blockSize, static_cast<std::align_val_t>(_alignment));
            cudaErrChk(cudaFree(_blocks[i].ptr));
#ifndef NDEBUG
            std::cout << "deallocated: " << _blocks[i].blockSize << "\t@: " << _blocks[i].ptr << '\n';
#endif
            _blocks[i].ptr = nullptr;
        }
    }
}

void* Pool::do_allocate(std::size_t bytes, std::size_t)
{
    // update max block size
    _maxBlockSize = std::max(bytes, _maxBlockSize);

    std::uint32_t i = _blocks.size();
    // search for free block
    if(!_freeBlocks.empty())
    {
        // there is at least a free block
        i = _freeBlocks.top();
        _freeBlocks.pop();
    }

    // check if we run out of blocks
    if(i == _blocks.size())
    {
        _blocks.push_back(details::block{});
    }

    // check if it was allocated before
    if(_blocks[i].ptr != nullptr)
    {
        // check if big enough
        if (bytes > _blocks[i].blockSize)
        {
            // realloc
            // ::operator delete(_blocks[i].ptr, _blocks[i].blockSize, static_cast<std::align_val_t>(_alignment));
            cudaErrChk(cudaFree(_blocks[i].ptr));
            // _blocks[i].ptr = ::operator new(_maxBlockSize, static_cast<std::align_val_t>(_alignment));
            cudaErrChk(cudaMalloc(reinterpret_cast<void**>(&_blocks[i].ptr), _maxBlockSize));
            _blocks[i].blockSize = _maxBlockSize;
#ifndef NDEBUG
            std::cout << "realloc block: " << i << "\tof size: " << _blocks[i].blockSize << '\n';
#endif
        }
        else
        {
#ifndef NDEBUG
            std::cout << "reuse block: " << i << "\tof size: " << _blocks[i].blockSize << '\n';
#endif
        }
        _blocks[i].isValid = true;
        return _blocks[i].ptr;
    }

    // if not allocate
    // auto ptr = ::operator new(_maxBlockSize, static_cast<std::align_val_t>(_alignment));
    void* ptr{nullptr};
    cudaErrChk(cudaMalloc(reinterpret_cast<void**>(&ptr), _maxBlockSize));
    _blocks[i].blockSize = _maxBlockSize;
#ifndef NDEBUG
    std::cout << "allocated: " << _blocks[i].blockSize  << "\t@: " << ptr << '\n';
#endif

    // store
    _blocks[i].ptr = ptr;
    _blocks[i].isValid = true;
    return ptr;
}

void Pool::do_deallocate(void* ptr, std::size_t bytes, std::size_t)
{
    assert(ptr != nullptr);

    for(std::uint32_t i = 0; i < _blocks.size(); ++i)
    {
        if(_blocks[i].ptr == ptr)
        {
            _blocks[i].isValid = false;
            _freeBlocks.push(i);
#ifndef NDEBUG
            std::cout << "freed block: " << i << "\tof size: " << _blocks[i].blockSize <<'\n';
#endif
            return;
        }
    }

}

bool Pool::do_is_equal(const QuICC::Memory::memory_resource& other) const noexcept
{
    if( this == &other ) return true;
    const auto* const op = dynamic_cast<const Pool*>(&other);
    return op != nullptr;
}

} // namespace Cuda
} // namespace Memory
} // namespace QuICC

