#include <iostream>
#include <complex>
#include <cassert>

#include "Graph/Shims/MlirShims.hpp"
#include "Graph/BackendsMap.hpp"
#include "Graph/Types.hpp"

using namespace QuICC::Graph;

/// @brief C Interface to MLIR for a transpose operator
/// cpu backend
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_201_complexf64_DCCSC3D_complexf64_DCCSC3D(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_201_complexf64_DCCSC3D_complexf64_DCCSC3D\n";
    #endif
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    assert(pIn->dataSize >= pOut->dataSize); // Input might be padded
    // op
    using namespace QuICC::Transpose::Cpu;
    using namespace QuICC::Transpose;
    using Tin = C_DCCSC3D_t;
    using Tout = C_DCCSC3D_t;
    using op_t = Op<Tout, Tin, p201_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> pointers[rank];
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> indices[rank];
    std::uint32_t lds = pIn->dataSize / pIn->cooSize;
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices, lds);
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a transpose operator
/// gpu backend
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_201_complexf64_DCCSC3DJIK_complexf64_DCCSC3D(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_201_complexf64_DCCSC3DJIK_complexf64_DCCSC3D\n";
    #endif
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    assert(pIn->dataSize == pOut->dataSize);
    assert(QuICC::Cuda::isDeviceMemory(pOut->data));
    assert(QuICC::Cuda::isDeviceMemory(pOut->pos));
    assert(QuICC::Cuda::isDeviceMemory(pOut->coo));
    assert(QuICC::Cuda::isDeviceMemory(pIn->data));
    assert(QuICC::Cuda::isDeviceMemory(pIn->pos));
    assert(QuICC::Cuda::isDeviceMemory(pIn->coo));
    // op
    using namespace QuICC::Transpose::Cuda;
    using namespace QuICC::Transpose;
    using Tin = C_DCCSC3D_t;
    using Tout = C_DCCSC3DJIK_t;
    using op_t = Op<Tout, Tin, p201_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> pointers[rank];
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> indices[rank];
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices);
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
};
#endif

/// @brief C Interface to MLIR for a transpose operator
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_120_complexf64_DCCSC3D_complexf64_DCCSC3D(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_120_complexf64_DCCSC3D_complexf64_DCCSC3D\n";
    #endif
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    assert(pIn->dataSize <= pOut->dataSize); // there could be padding for FFT
    // op
    using namespace QuICC::Transpose::Cpu;
    using namespace QuICC::Transpose;
    using Tin = C_DCCSC3D_t;
    using Tout = C_DCCSC3D_t;
    using op_t = Op<Tout, Tin, p120_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> pointers[rank];
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> indices[rank];
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices);
    std::uint32_t lds = pOut->dataSize / pOut->cooSize;
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices, lds);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a transpose operator
/// gpu backend
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_120_complexf64_DCCSC3D_complexf64_DCCSC3DJIK(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_120_complexf64_DCCSC3D_complexf64_DCCSC3DJIK\n";
    #endif
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    assert(pIn->dataSize == pOut->dataSize);
    assert(QuICC::Cuda::isDeviceMemory(pOut->data));
    assert(QuICC::Cuda::isDeviceMemory(pOut->pos));
    assert(QuICC::Cuda::isDeviceMemory(pOut->coo));
    assert(QuICC::Cuda::isDeviceMemory(pIn->data));
    assert(QuICC::Cuda::isDeviceMemory(pIn->pos));
    assert(QuICC::Cuda::isDeviceMemory(pIn->coo));
    // op
    using namespace QuICC::Transpose::Cuda;
    using namespace QuICC::Transpose;
    using Tout = C_DCCSC3D_t;
    using Tin = C_DCCSC3DJIK_t;
    using op_t = Op<Tout, Tin, p120_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> pointers[rank];
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> indices[rank];
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices);
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
};
#endif

/// @brief C Interface to MLIR for a transpose operator
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_201_complexf64_DCCSC3D_complexf64_S1CLCSC3D(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_201_complexf64_DCCSC3D_complexf64_S1CLCSC3D\n";
    #endif
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    assert(pIn->dataSize == pOut->dataSize);
    // Op
    using namespace QuICC::Transpose::Cpu;
    using namespace QuICC::Transpose;
    using Tin = C_S1CLCSC3D_t;
    using Tout = C_DCCSC3D_t;
    using op_t = Op<Tout, Tin, p201_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> pointers[rank];
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> indices[rank];
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices);
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a transpose operator
/// gpu backend
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_201_complexf64_DCCSC3DJIK_complexf64_S1CLCSC3DJIK(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_201_complexf64_DCCSC3DJIK_complexf64_S1CLCSC3DJIK\n";
    #endif
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    assert(pIn->dataSize == pOut->dataSize);
    assert(QuICC::Cuda::isDeviceMemory(pOut->data));
    assert(QuICC::Cuda::isDeviceMemory(pOut->pos));
    assert(QuICC::Cuda::isDeviceMemory(pOut->coo));
    assert(QuICC::Cuda::isDeviceMemory(pIn->data));
    assert(QuICC::Cuda::isDeviceMemory(pIn->pos));
    assert(QuICC::Cuda::isDeviceMemory(pIn->coo));
    // Op
    using namespace QuICC::Transpose::Cuda;
    using namespace QuICC::Transpose;
    using Tin = C_S1CLCSC3DJIK_t;
    using Tout = C_DCCSC3DJIK_t;
    using op_t = Op<Tout, Tin, p201_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> pointers[rank];
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> indices[rank];
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices);
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
};
#endif

/// @brief C Interface to MLIR for a transpose operator
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_120_complexf64_S1CLCSC3D_complexf64_DCCSC3D(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_120_complexf64_S1CLCSC3D_complexf64_DCCSC3D\n";
    #endif
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    assert(pIn->dataSize == pOut->dataSize);
    // Op
    using namespace QuICC::Transpose::Cpu;
    using namespace QuICC::Transpose;
    using Tin = C_DCCSC3D_t;
    using Tout = C_S1CLCSC3D_t;
    using op_t = Op<Tout, Tin, p120_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> pointers[rank];
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> indices[rank];
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices);
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a transpose operator
/// gpu backend
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_120_complexf64_S1CLCSC3DJIK_complexf64_DCCSC3DJIK(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_120_complexf64_S1CLCSC3DJIK_complexf64_DCCSC3DJIK\n";
    #endif
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    assert(pIn->dataSize == pOut->dataSize);
    assert(QuICC::Cuda::isDeviceMemory(pOut->data));
    assert(QuICC::Cuda::isDeviceMemory(pOut->pos));
    assert(QuICC::Cuda::isDeviceMemory(pOut->coo));
    assert(QuICC::Cuda::isDeviceMemory(pIn->data));
    assert(QuICC::Cuda::isDeviceMemory(pIn->pos));
    assert(QuICC::Cuda::isDeviceMemory(pIn->coo));
    // Op
    using namespace QuICC::Transpose::Cuda;
    using namespace QuICC::Transpose;
    using Tin = C_DCCSC3DJIK_t;
    using Tout = C_S1CLCSC3DJIK_t;
    using op_t = Op<Tout, Tin, p120_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> pointers[rank];
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> indices[rank];
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices);
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
};
#endif
