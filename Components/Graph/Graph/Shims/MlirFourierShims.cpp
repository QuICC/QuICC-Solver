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
extern "C" inline void _ciface_quiccir_fr_prj_f64_DCCSC3D_complexf64_DCCSC3D_cpu(void* obj, view3_t* pUval, view3_cd_t* pUmod)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_fr_prj_f64_DCCSC3D_complexf64_DCCSC3D_cpu\n";
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
extern "C" inline void _ciface_quiccir_fr_prj_f64_DCCSC3D_complexf64_DCCSC3D_gpu(void* obj, view3_t* pUval, view3_cd_t* pUmod)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_fr_prj_f64_DCCSC3D_complexf64_DCCSC3D_gpu\n";
    #endif
    assert(obj != nullptr);
    assert(pUval != nullptr);
    assert(pUmod != nullptr);
    assert(QuICC::Cuda::isDeviceMemory(pUmod->data));
    assert(QuICC::Cuda::isDeviceMemory(pUmod->pos));
    assert(QuICC::Cuda::isDeviceMemory(pUmod->coo));
    assert(QuICC::Cuda::isDeviceMemory(pUval->data));
    assert(QuICC::Cuda::isDeviceMemory(pUval->pos));
    assert(QuICC::Cuda::isDeviceMemory(pUval->coo));
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
extern "C" void _ciface_quiccir_fr_prj_f64_DCCSC3D_complexf64_DCCSC3D(void* obj, view3_t* pUval, view3_cd_t* pUmod)
{
    #ifdef QUICC_HAS_CUDA_BACKEND
    assert(pUval != nullptr);
    if (!QuICC::Cuda::isDeviceMemory(pUval->data))
    {
        assert(!QuICC::Cuda::isDeviceMemory(pUmod->data));
        _ciface_quiccir_fr_prj_f64_DCCSC3D_complexf64_DCCSC3D_cpu(obj, pUval, pUmod);
    }
    else
    {
        _ciface_quiccir_fr_prj_f64_DCCSC3D_complexf64_DCCSC3D_gpu(obj, pUval, pUmod);
    }
    #else
    _ciface_quiccir_fr_prj_f64_DCCSC3D_complexf64_DCCSC3D_cpu(obj, pUval, pUmod);
    #endif
}

/// @brief C Interface to MLIR for a fr int operator
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_fr_int_complexf64_DCCSC3D_f64_DCCSC3D_cpu(void* obj, view3_cd_t* pUmod, view3_t* pUval)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_fr_int_complexf64_DCCSC3D_f64_DCCSC3D_cpu\n";
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
extern "C" void _ciface_quiccir_fr_int_complexf64_DCCSC3D_f64_DCCSC3D_gpu(void* obj, view3_cd_t* pUmod, view3_t* pUval)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_fr_int_complexf64_DCCSC3D_f64_DCCSC3D_gpu\n";
    #endif
    assert(obj != nullptr);
    assert(pUval != nullptr);
    assert(pUmod != nullptr);
    assert(QuICC::Cuda::isDeviceMemory(pUmod->data));
    assert(QuICC::Cuda::isDeviceMemory(pUmod->pos));
    assert(QuICC::Cuda::isDeviceMemory(pUmod->coo));
    assert(QuICC::Cuda::isDeviceMemory(pUval->data));
    assert(QuICC::Cuda::isDeviceMemory(pUval->pos));
    assert(QuICC::Cuda::isDeviceMemory(pUval->coo));
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
extern "C" void _ciface_quiccir_fr_int_complexf64_DCCSC3D_f64_DCCSC3D(void* obj, view3_cd_t* pUmod, view3_t* pUval)
{
    #ifdef QUICC_HAS_CUDA_BACKEND
    assert(pUval != nullptr);
    if (!QuICC::Cuda::isDeviceMemory(pUval->data))
    {
        assert(!QuICC::Cuda::isDeviceMemory(pUmod->data));
        _ciface_quiccir_fr_int_complexf64_DCCSC3D_f64_DCCSC3D_cpu(obj, pUmod, pUval);
    }
    else
    {
        _ciface_quiccir_fr_int_complexf64_DCCSC3D_f64_DCCSC3D_gpu(obj, pUmod, pUval);
    }
    #else
    _ciface_quiccir_fr_int_complexf64_DCCSC3D_f64_DCCSC3D_cpu(obj, pUmod, pUval);
    #endif
}
