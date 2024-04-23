#include <iostream>
#include <complex>
#include <cassert>

#include "Graph/MlirShims.hpp"
#include "Graph/BackendsMap.hpp"
#include "View/View.hpp"

using namespace QuICC::Graph;

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_R_DCCSC3D_t(view3_t* pBuffer)
{
    // Reset Meta
    pBuffer->coo = nullptr;
    pBuffer->cooSize = 0;
    pBuffer->pos = nullptr;
    pBuffer->posSize = 0;
    // Dealloc
    assert(pBuffer->data != nullptr);
    std::size_t sizeByte = sizeof(double) * pBuffer->dataSize;
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
    std::cout << "_ciface_quiccir_dealloc_R_DCCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_C_DCCSC3D_t(view3_t* pBuffer)
{
    // Reset Meta
    pBuffer->coo = nullptr;
    pBuffer->cooSize = 0;
    pBuffer->pos = nullptr;
    pBuffer->posSize = 0;
    // Dealloc
    assert(pBuffer->data != nullptr);
    std::size_t sizeByte = sizeof(std::complex<double>) * pBuffer->dataSize;
    ::operator delete(pBuffer->data, sizeByte, static_cast<std::align_val_t>(sizeof(std::complex<double>)));
    pBuffer->dataSize = 0;
    pBuffer->data = nullptr;
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_C_DCCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_C_S1CLCSC3D_t(view3_t* pBuffer)
{
    // Reset Meta
    pBuffer->coo = nullptr;
    pBuffer->cooSize = 0;
    pBuffer->pos = nullptr;
    pBuffer->posSize = 0;
    // Dealloc
    assert(pBuffer->data != nullptr);
    std::size_t sizeByte = sizeof(std::complex<double>) * pBuffer->dataSize;
    ::operator delete(pBuffer->data, sizeByte, static_cast<std::align_val_t>(sizeof(std::complex<double>)));
    pBuffer->dataSize = 0;
    pBuffer->data = nullptr;
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_C_S1CLCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};
