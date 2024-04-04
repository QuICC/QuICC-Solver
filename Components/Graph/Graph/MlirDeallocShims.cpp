#include <iostream>
#include <complex>
#include <cassert>
#include <Quiccir-c/Utils.h>

#include "Graph/MlirShims.hpp"
#include "Graph/BackendsMap.hpp"
#include "View/View.hpp"

using namespace QuICC::Graph;

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_R_DCCSC3D_t(view3_t* pBuffer)
{
    // meta
    pBuffer->coo = nullptr;
    pBuffer->cooSize = 0;
    pBuffer->pos = nullptr;
    pBuffer->posSize = 0;
    // dealloc
    assert(pBuffer->data != nullptr);
    std::size_t sizeByte = sizeof(double) * pBuffer->dataSize;
    ::operator delete(pBuffer->data, sizeByte, static_cast<std::align_val_t>(sizeof(double)));
    pBuffer->dataSize = 0;
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_R_DCCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_C_DCCSC3D_t(view3_t* pBuffer)
{
    // meta
    pBuffer->coo = nullptr;
    pBuffer->cooSize = 0;
    pBuffer->pos = nullptr;
    pBuffer->posSize = 0;
    // dealloc
    assert(pBuffer->data != nullptr);
    std::size_t sizeByte = sizeof(std::complex<double>) * pBuffer->dataSize;
    ::operator delete(pBuffer->data, sizeByte, static_cast<std::align_val_t>(sizeof(std::complex<double>)));
    pBuffer->dataSize = 0;
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_C_DCCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_C_S1CLCSC3D_t(view3_t* pBuffer)
{
    // meta
    pBuffer->coo = nullptr;
    pBuffer->cooSize = 0;
    pBuffer->pos = nullptr;
    pBuffer->posSize = 0;
    // dealloc
    assert(pBuffer->data != nullptr);
    std::size_t sizeByte = sizeof(std::complex<double>) * pBuffer->dataSize;
    ::operator delete(pBuffer->data, sizeByte, static_cast<std::align_val_t>(sizeof(std::complex<double>)));
    pBuffer->dataSize = 0;
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_C_S1CLCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};
