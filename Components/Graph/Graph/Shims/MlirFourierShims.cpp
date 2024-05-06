#include <iostream>
#include <complex>
#include <cassert>

#include "Graph/BackendsMap.hpp"
#include "Graph/Shims/MlirShims.hpp"
#include "Graph/Types.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "Cuda/CudaUtil.hpp"
#endif

using namespace QuICC::Graph;

/// @brief C Interface to MLIR for a fr prj operator, cpu backend
/// @param op
/// @param uval
/// @param umod
extern "C" inline void _ciface_quiccir_fr_prj_R_DCCSC3D_t_C_DCCSC3D_t_cpu(void* obj, view3_t* pUval, view3_cd_t* pUmod)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_fr_prj_R_DCCSC3D_t_C_DCCSC3D_t_cpu\n";
    #endif
    assert(obj != nullptr);
    assert(pUval != nullptr);
    assert(pUmod != nullptr);
    // op
    using namespace QuICC::Transform::Fourier;
    using backend_t = QuICC::Graph::viewCpu_t;
    using Tin = C_DCCSC3D_t;
    using Tout = R_DCCSC3D_t;
    using backendFft_t = Fft_t<backend_t, Tout, Tin>;
    using backendDiff_t = MixedDiff_t<backend_t, Tin, 0, bwd_t,
        QuICC::Transform::Fourier::none_m>;
    using op_t = Mixed::Projector::DOp<Tout, Tin, backendFft_t,   backendDiff_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    ViewBase<std::uint32_t> pointers[rank];
    pointers[1] = ViewBase<std::uint32_t>(pUmod->pos, pUmod->posSize);
    ViewBase<std::uint32_t> indices[rank];
    indices[1] = ViewBase<std::uint32_t>(pUmod->coo, pUmod->cooSize);
    Tin viewMod(pUmod->data, pUmod->dataSize, pUmod->dims, pointers, indices);
    Tout viewVal(pUval->data, pUval->dataSize, pUval->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewVal, viewMod);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a fr prj operator, cuda backend
/// @param op
/// @param uval
/// @param umod
extern "C" inline void _ciface_quiccir_fr_prj_R_DCCSC3D_t_C_DCCSC3D_t_gpu(void* obj, view3_t* pUval, view3_cd_t* pUmod)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_fr_prj_R_DCCSC3D_t_C_DCCSC3D_t_gpu\n";
    #endif
    assert(obj != nullptr);
    assert(pUval != nullptr);
    assert(pUmod != nullptr);
    // op
    using namespace QuICC::Transform::Fourier;
    using backend_t = QuICC::Graph::viewGpu_t;
    using Tin = C_DCCSC3D_t;
    using Tout = R_DCCSC3D_t;
    using backendFft_t = Fft_t<backend_t, Tout, Tin>;
    using backendDiff_t = MixedDiff_t<backend_t, Tin, 0, bwd_t,
        QuICC::Transform::Fourier::none_m>;
    using op_t = Mixed::Projector::DOp<Tout, Tin, backendFft_t,   backendDiff_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    ViewBase<std::uint32_t> pointers[rank];
    pointers[1] = ViewBase<std::uint32_t>(pUmod->pos, pUmod->posSize);
    ViewBase<std::uint32_t> indices[rank];
    indices[1] = ViewBase<std::uint32_t>(pUmod->coo, pUmod->cooSize);
    Tin viewMod(pUmod->data, pUmod->dataSize, pUmod->dims, pointers, indices);
    Tout viewVal(pUval->data, pUval->dataSize, pUval->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewVal, viewMod);
};
#endif

/// @brief C Interface to MLIR for a fr prj operator
/// @param op
/// @param uval
/// @param umod
extern "C" void _ciface_quiccir_fr_prj_R_DCCSC3D_t_C_DCCSC3D_t(void* obj, view3_t* pUval, view3_cd_t* pUmod)
{
    #ifdef QUICC_HAS_CUDA_BACKEND
    assert(pUval != nullptr);
    if (!QuICC::Cuda::isDeviceMemory(pUval->data))
    {
        assert(!QuICC::Cuda::isDeviceMemory(pUmod->data));
        _ciface_quiccir_fr_prj_R_DCCSC3D_t_C_DCCSC3D_t_cpu(obj, pUval, pUmod);
    }
    else
    {
        _ciface_quiccir_fr_prj_R_DCCSC3D_t_C_DCCSC3D_t_gpu(obj, pUval, pUmod);
    }
    #else
    _ciface_quiccir_fr_prj_R_DCCSC3D_t_C_DCCSC3D_t_cpu(obj, pUval, pUmod);
    #endif
}

/// @brief C Interface to MLIR for a fr int operator
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_fr_int_C_DCCSC3D_t_R_DCCSC3D_t_cpu(void* obj, view3_cd_t* pUmod, view3_t* pUval)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_fr_int_C_DCCSC3D_t_R_DCCSC3D_t_cpu\n";
    #endif
    assert(obj != nullptr);
    assert(pUval != nullptr);
    assert(pUmod != nullptr);
    // op
    using namespace QuICC::Transform::Fourier;
    using backend_t = QuICC::Graph::viewCpu_t;
    using Tin = R_DCCSC3D_t;
    using Tout = C_DCCSC3D_t;
    using backendFft_t = Fft_t<backend_t, Tout, Tin>;
    using backendDiff_t = MixedDiff_t<backend_t, Tout, 0, fwd_t,
        QuICC::Transform::Fourier::none_m>;
    using op_t = Mixed::Integrator::DOp<Tout, Tin, backendFft_t,   backendDiff_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    ViewBase<std::uint32_t> pointers[rank];
    pointers[1] = ViewBase<std::uint32_t>(pUmod->pos, pUmod->posSize);
    ViewBase<std::uint32_t> indices[rank];
    indices[1] = ViewBase<std::uint32_t>(pUmod->coo, pUmod->cooSize);
    Tout viewMod(pUmod->data, pUmod->dataSize, pUmod->dims, pointers, indices);
    Tin viewVal(pUval->data, pUval->dataSize, pUval->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewMod, viewVal);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a fr int operator
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_fr_int_C_DCCSC3D_t_R_DCCSC3D_t_gpu(void* obj, view3_cd_t* pUmod, view3_t* pUval)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_fr_int_C_DCCSC3D_t_R_DCCSC3D_t_gpu\n";
    #endif
    assert(obj != nullptr);
    assert(pUval != nullptr);
    assert(pUmod != nullptr);
    // op
    using namespace QuICC::Transform::Fourier;
    using backend_t = QuICC::Graph::viewGpu_t;
    using Tin = R_DCCSC3D_t;
    using Tout = C_DCCSC3D_t;
    using backendFft_t = Fft_t<backend_t, Tout, Tin>;
    using backendDiff_t = MixedDiff_t<backend_t, Tout, 0, fwd_t,
        QuICC::Transform::Fourier::none_m>;
    using op_t = Mixed::Integrator::DOp<Tout, Tin, backendFft_t,   backendDiff_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    ViewBase<std::uint32_t> pointers[rank];
    pointers[1] = ViewBase<std::uint32_t>(pUmod->pos, pUmod->posSize);
    ViewBase<std::uint32_t> indices[rank];
    indices[1] = ViewBase<std::uint32_t>(pUmod->coo, pUmod->cooSize);
    Tout viewMod(pUmod->data, pUmod->dataSize, pUmod->dims, pointers, indices);
    Tin viewVal(pUval->data, pUval->dataSize, pUval->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewMod, viewVal);
};
#endif

/// @brief C Interface to MLIR for a fr int operator
/// @param op
/// @param uval
/// @param umod
extern "C" void _ciface_quiccir_fr_int_C_DCCSC3D_t_R_DCCSC3D_t(void* obj, view3_cd_t* pUmod, view3_t* pUval)
{
    #ifdef QUICC_HAS_CUDA_BACKEND
    assert(pUval != nullptr);
    if (!QuICC::Cuda::isDeviceMemory(pUval->data))
    {
        assert(!QuICC::Cuda::isDeviceMemory(pUmod->data));
        _ciface_quiccir_fr_int_C_DCCSC3D_t_R_DCCSC3D_t_cpu(obj, pUmod, pUval);
    }
    else
    {
        _ciface_quiccir_fr_int_C_DCCSC3D_t_R_DCCSC3D_t_gpu(obj, pUmod, pUval);
    }
    #else
    _ciface_quiccir_fr_int_C_DCCSC3D_t_R_DCCSC3D_t_cpu(obj, pUmod, pUval);
    #endif
}

/// @brief C Interface to MLIR for an allocator
/// The name fully describes how the buffer needs to be allocated
/// since it includes the producer op and buffer.
/// @param pNewBuffer buffer to be allocated
/// @param pProdBuffer producer buffer
extern "C" void _ciface_quiccir_alloc_fr_prj_R_DCCSC3D_t_C_DCCSC3D_t(view3_t* pNewBuffer, view3_cd_t* pProdBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_fr_prj_R_DCCSC3D_t_C_DCCSC3D_t\n";
    #endif
    // Slice (phys = mods*2-1)
    // Projector, NewBuffer is in physical space
    assert(pProdBuffer->dims[0] == std::floor(pNewBuffer->dims[0]/2)+1);
    assert(pNewBuffer->dims[1] == pProdBuffer->dims[1]);
    // Layers
    assert(pNewBuffer->dims[2] == pProdBuffer->dims[2]);
    // Reuse meta
    assert(pProdBuffer->pos != nullptr);
    assert(pProdBuffer->coo != nullptr);
    pNewBuffer->pos = pProdBuffer->pos;
    pNewBuffer->posSize = pProdBuffer->posSize;
    pNewBuffer->coo = pProdBuffer->coo;
    pNewBuffer->cooSize = pProdBuffer->cooSize;
    // Alloc buffer
    pNewBuffer->dataSize = pNewBuffer->dims[0] * pNewBuffer->cooSize;
    details::alloc_ptr(&pNewBuffer->data, pNewBuffer->dataSize, pProdBuffer->data);
};

extern "C" void _ciface_quiccir_alloc_fr_int_C_DCCSC3D_t_R_DCCSC3D_t(view3_cd_t* pNewBuffer, view3_t* pProdBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_fr_int_C_DCCSC3D_t_R_DCCSC3D_t\n";
    #endif
    // Slice (phys = mods*2-1)
    // Integrator, NewBuffer is in modal space
    assert(pNewBuffer->dims[0] == std::floor(pProdBuffer->dims[0]/2)+1);
    assert(pNewBuffer->dims[1] == pProdBuffer->dims[1]);
    // Layers
    assert(pNewBuffer->dims[2] == pProdBuffer->dims[2]);
    // Reuse meta
    assert(pProdBuffer->pos != nullptr);
    assert(pProdBuffer->coo != nullptr);
    pNewBuffer->pos = pProdBuffer->pos;
    pNewBuffer->posSize = pProdBuffer->posSize;
    pNewBuffer->coo = pProdBuffer->coo;
    pNewBuffer->cooSize = pProdBuffer->cooSize;
    // Alloc buffer
    pNewBuffer->dataSize = pNewBuffer->dims[0] * pNewBuffer->cooSize;
    details::alloc_ptr(&pNewBuffer->data, pNewBuffer->dataSize, pProdBuffer->data);
};
