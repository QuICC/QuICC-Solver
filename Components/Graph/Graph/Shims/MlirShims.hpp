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

    /// @brief compute pointers and index meta data for fully populated
    /// tensor as output of a Transpose operation
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
    /// tensor as ouput of a Transpose operation between JW and AL
    /// @tparam C_S1CLCSC3D_t new buffer type
    /// @tparam C_DCCSC3D_t producer type
    /// @param dims output dimensions
    /// @return
    template<>
    ptrAndIdx denseTransposePtrAndIdx<C_S1CLCSC3D_t, C_DCCSC3D_t>(const std::array<std::uint32_t, 3> dims);

    /// @brief compute pointers and index meta data for fully populated
    /// tensor as ouput of a Transpose operation between AL and JW
    /// @tparam C_DCCSC3DJIK_t new buffer type
    /// @tparam C_S1CLCSC3DJIK_t producer type
    /// @param dims output dimensions
    /// @return
    template<>
    ptrAndIdx denseTransposePtrAndIdx<C_DCCSC3DJIK_t, C_S1CLCSC3DJIK_t>(const std::array<std::uint32_t, 3> dims);

    /// @brief compute pointers and index meta data for fully populated
    /// tensor as ouput of a Transpose operation between JW and AL
    /// @tparam C_S1CLCSC3DJIK_t new buffer type
    /// @tparam C_DCCSC3DJIK_t producer type
    /// @param dims output dimensions
    /// @return
    template<>
    ptrAndIdx denseTransposePtrAndIdx<C_S1CLCSC3DJIK_t, C_DCCSC3DJIK_t>(const std::array<std::uint32_t, 3> dims);

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
        auto alignment = static_cast<std::align_val_t>(alignof(std::max_align_t));
        *newPtr = reinterpret_cast<Tnew*>(::operator new(sizeByte, alignment));
    }
    #ifdef QUICC_HAS_CUDA_BACKEND
    else
    {
        cudaErrChk(cudaMalloc(reinterpret_cast<void**>(newPtr), sizeByte));
    }
    #endif
    #ifndef NDEBUG
    std::cout << "alloc, bytes: " << sizeByte << "\tis cpu: " << isCpuMem << "\tptr: " << *newPtr <<'\n';
    #endif
}


/// @brief generic view descriptor deallocator
/// @tparam Tdata scalar type of view descriptor
/// @param pBuffer
template<class Tdata>
void dealloc_viewDescriptor(ViewDescriptor<Tdata, std::uint32_t, 3>* pBuffer)
{
    // Reset Meta
    pBuffer->coo = nullptr;
    pBuffer->cooSize = 0;
    pBuffer->pos = nullptr;
    pBuffer->posSize = 0;
    // Dealloc
    assert(pBuffer->data != nullptr);
    std::size_t sizeByte = sizeof(Tdata) * pBuffer->dataSize;
    // Check memory space
    bool isCpuMem = true;
    #ifndef NDEBUG
    std::cout << "dealloc, bytes: " << sizeByte << "\tis cpu: " << isCpuMem << "\tptr: " << pBuffer->data <<'\n';
    #endif
    #ifdef QUICC_HAS_CUDA_BACKEND
    isCpuMem = !QuICC::Cuda::isDeviceMemory(pBuffer->data);
    #endif
    if (isCpuMem)
    {
        auto alignment = static_cast<std::align_val_t>(alignof(std::max_align_t));
        ::operator delete(pBuffer->data, sizeByte, alignment);
    }
    #ifdef QUICC_HAS_CUDA_BACKEND
    else
    {
        cudaErrChk(cudaFree(pBuffer->data));
    }
    #endif
    pBuffer->dataSize = 0;
    pBuffer->data = nullptr;
}

#ifdef QUICC_HAS_CUDA_BACKEND
std::uint32_t getSizeS1CLCSC3DJIK(const std::uint32_t* ptr, const std::uint32_t size, const std::uint32_t lds);
#endif

} // namespace details
} // namespace Graph
} // namespace QuICC
