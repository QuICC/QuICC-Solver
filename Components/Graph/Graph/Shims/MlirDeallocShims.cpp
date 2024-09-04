#include <iostream>
#include <complex>
#include <cassert>

#include "Graph/Shims/MlirShims.hpp"
#include "View/View.hpp"

using namespace QuICC::Graph;

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_f64_DCCSC3D(view3_t* pBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_f64_DCCSC3D\n";
    #endif
    details::dealloc_viewDescriptor(pBuffer);
};

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_complexf64_DCCSC3D(view3_cd_t* pBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_complexf64_DCCSC3D\n";
    #endif
    details::dealloc_viewDescriptor(pBuffer);
};

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_complexf64_DCCSC3DJIK(view3_cd_t* pBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_complexf64_DCCSC3DJIK\n";
    #endif
    details::dealloc_viewDescriptor(pBuffer);
};

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_complexf64_S1CLCSC3D(view3_cd_t* pBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_complexf64_S1CLCSC3D\n";
    #endif
    details::dealloc_viewDescriptor(pBuffer);
};

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_complexf64_S1CLCSC3DJIK(view3_cd_t* pBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_complexf64_S1CLCSC3DJIK\n";
    #endif
    details::dealloc_viewDescriptor(pBuffer);
};
