#include <iostream>
#include <complex>
#include <cassert>

#include "Graph/Shims/MlirShims.hpp"
#include "Graph/BackendsMap.hpp"
#include "View/View.hpp"

using namespace QuICC::Graph;

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_R_DCCSC3D_t(view3_t* pBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_R_DCCSC3D_t\n";
    #endif
    details::dealloc_viewDescriptor(pBuffer);
};

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_C_DCCSC3D_t(view3_cd_t* pBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_C_DCCSC3D_t\n";
    #endif
    details::dealloc_viewDescriptor(pBuffer);
};

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_C_DCCSC3DJIK_t(view3_cd_t* pBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_C_DCCSC3DJIK_t\n";
    #endif
    details::dealloc_viewDescriptor(pBuffer);
};

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_C_S1CLCSC3D_t(view3_cd_t* pBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_C_S1CLCSC3D_t\n";
    #endif
    details::dealloc_viewDescriptor(pBuffer);
};

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_C_S1CLCSC3DJIK_t(view3_cd_t* pBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_C_S1CLCSC3DJIK_t\n";
    #endif
    details::dealloc_viewDescriptor(pBuffer);
};
