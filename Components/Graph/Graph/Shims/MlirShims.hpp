/**
 * @file MlirShims.hpp
 * @brief
 */
#pragma once

// External includes
//
#include <complex>
#include <unordered_map>
#include <iostream>

// Project includes
//
#include "Graph/Types.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"

namespace QuICC
{
namespace Graph
{
    template<typename Tdata, typename Tmeta, std::size_t N>
    struct ViewDescriptor {
    Tmeta dims[N];
    Tmeta *pos;
    Tmeta posSize;
    Tmeta *coo;
    Tmeta cooSize;
    Tdata *data;
    Tmeta dataSize;
    };

    using view3_t = ViewDescriptor<double, std::uint32_t, 3>;
    using view3_cd_t = ViewDescriptor<std::complex<double>, std::uint32_t, 3>;


    struct ptrAndIdx
    {
        std::vector<std::uint32_t> ptr;
        std::vector<std::uint32_t> idx;
    };

    template <class TagTo, class TagFrom>
    ptrAndIdx denseTransposePtrAndIdx(const std::array<std::uint32_t, 3> dims)
    {
        assert(false && "not implemented");
        return ptrAndIdx{};
    }

    /// @brief compute pointers and index meta data for fully populated
    /// tensor as ouput of a Transpose operation between AL and JW
    /// @tparam C_DCCSC3D_t new buffer type
    /// @tparam C_S1CLCSC3D_t producer type
    /// @param dims output dimensions
    /// @return
    template<>
    ptrAndIdx denseTransposePtrAndIdx<C_DCCSC3D_t, C_S1CLCSC3D_t>(const std::array<std::uint32_t, 3> dims);

    /// @brief compute pointers and index meta data for fully populated
    /// tensor as ouput of a Transpose operation between AL and JW
    /// @tparam C_S1CLCSC3D_t new buffer type
    /// @tparam C_DCCSC3D_t producer type
    /// @param dims output dimensions
    /// @return
    template<>
    ptrAndIdx denseTransposePtrAndIdx<C_S1CLCSC3D_t, C_DCCSC3D_t>(const std::array<std::uint32_t, 3> dims);

namespace details
{

/// @brief generic view descriptor allocator
/// @tparam T view descriptor type
/// @param pBuffer
template<class Tnew, class Tprod>
void alloc_ptr(Tnew** newPtr, const size_t size, const Tprod* prodPtr)
{
    std::size_t sizeByte = sizeof(Tnew) * size;
    // Check memory space
    bool isCpuMem = true;
    #ifdef QUICC_HAS_CUDA_BACKEND
    isCpuMem = !QuICC::Cuda::isDeviceMemory(prodPtr);
    #endif
    if (isCpuMem)
    {
        *newPtr = reinterpret_cast<Tnew*>(::operator new(sizeByte, static_cast<std::align_val_t>(sizeof(Tnew))));
    }
    #ifdef QUICC_HAS_CUDA_BACKEND
    else
    {
        cudaErrChk(cudaMalloc(reinterpret_cast<void**>(newPtr), sizeByte));
    }
    #endif
    #ifndef NDEBUG
    std::cout << "alloc, bytes: " << sizeByte << '\t' << isCpuMem <<'\n';
    #endif
}


/// @brief generic view descriptor deallocator
/// @tparam T view descriptor type
/// @param pBuffer
template<class T>
void dealloc_viewDescriptor(T* pBuffer)
{
    // Reset Meta
    pBuffer->coo = nullptr;
    pBuffer->cooSize = 0;
    pBuffer->pos = nullptr;
    pBuffer->posSize = 0;
    // Dealloc
    assert(pBuffer->data != nullptr);
    std::size_t sizeByte = sizeof(decltype(*pBuffer->data)) * pBuffer->dataSize;
    // Check memory space
    bool isCpuMem = true;
    #ifdef QUICC_HAS_CUDA_BACKEND
    isCpuMem = !QuICC::Cuda::isDeviceMemory(pBuffer->data);
    #endif
    if (isCpuMem)
    {
        ::operator delete(pBuffer->data, sizeByte, static_cast<std::align_val_t>(sizeof(double)));
    }
    #ifdef QUICC_HAS_CUDA_BACKEND
    else
    {
        cudaErrChk(cudaFree(pBuffer->data));
    }
    #endif
    pBuffer->dataSize = 0;
    pBuffer->data = nullptr;
    #ifndef NDEBUG
    std::cout << "dealloc, bytes: " << sizeByte << '\t' << isCpuMem <<'\n';
    #endif
}

} // namespace details
} // namespace Graph
} // namespace QuICC


