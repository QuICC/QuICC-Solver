#include <iostream>
#include <complex>
#include <cassert>

#include "Graph/Shims/MlirShims.hpp"
#include "ViewOps/Slicewise/Op.hpp"
#include "ViewOps/Slicewise/Functors.hpp"
#include "Graph/Types.hpp"


using namespace QuICC::Graph;

/// @brief C Interface to MLIR for a mul const operator
/// @param obj
/// @param pRet
/// @param pU0
/// @param pU1
/// @param pU2
/// @param pV0
/// @param pV1
/// @param pV2
extern "C" void _ciface_quiccir_mul_const_buoyancy_f64_DCCSC3D_f64_DCCSC3D_cpu(void* obj,
    ViewDescriptor<double, std::uint32_t, 3>* pRet,
    ViewDescriptor<double, std::uint32_t, 3>* pRhs)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_mul_const_buoyancy_f64_DCCSC3D_f64_DCCSC3D_cpu\n";
    #endif
    assert(obj != nullptr);
    assert(pRet != nullptr);
    assert(pRhs != nullptr);
    // op
    using namespace QuICC::Slicewise::Cpu;
    using namespace QuICC::Slicewise;
    using T = QuICC::View::View<double, QuICC::View::DCCSC3D>;
    using op_t = Op<MulRFunctor<double>, T, T>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    ViewBase<std::uint32_t> pointers[rank];
    pointers[1] = ViewBase<std::uint32_t>(pRhs->pos, pRhs->posSize);
    ViewBase<std::uint32_t> indices[rank];
    indices[1] = ViewBase<std::uint32_t>(pRhs->coo, pRhs->cooSize);
    T viewRhs(pRhs->data, pRhs->dataSize, pRhs->dims, pointers, indices);
    T viewRet(pRet->data, pRet->dataSize, pRet->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewRet, viewRhs);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a binary mul const operator
/// gpu backend
/// @param obj
/// @param pRet
/// @param pRhs
extern "C" void _ciface_quiccir_mul_const_buoyancy_f64_DCCSC3D_f64_DCCSC3D_gpu(void* obj,
    ViewDescriptor<double, std::uint32_t, 3>* pRet,
    ViewDescriptor<double, std::uint32_t, 3>* pRhs)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_mul_const_buoyancy_f64_DCCSC3D_f64_DCCSC3D_gpu\n";
    #endif
    assert(obj != nullptr);
    assert(pRet != nullptr);
    assert(pRhs != nullptr);
    // op
    using namespace QuICC::Pointwise::Cuda;
    using namespace QuICC::Pointwise;
    using T = QuICC::View::View<double, QuICC::View::DCCSC3D>;
    using op_t = Op<MulRFunctor<double>, T, T>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    ViewBase<std::uint32_t> pointers[rank];
    pointers[1] = ViewBase<std::uint32_t>(pLhs->pos, pLhs->posSize);
    ViewBase<std::uint32_t> indices[rank];
    indices[1] = ViewBase<std::uint32_t>(pLhs->coo, pLhs->cooSize);
    T viewRhs(reinterpret_cast<cuda::double*>(pRhs->data), pRhs->dataSize, pRhs->dims, pointers, indices);
    T viewRet(reinterpret_cast<cuda::double*>(pRet->data), pRet->dataSize, pRet->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewRet, viewRhs);
};
#endif

/// @brief C Interface to MLIR for a mul const operator
/// @param op
/// @param uval
/// @param umod
extern "C" void _ciface_quiccir_mul_const_buoyancy_f64_DCCSC3D_f64_DCCSC3D(void* obj, view3_cd_t* pUmod, view3_t* pUval)
{
    #ifdef QUICC_HAS_CUDA_BACKEND
    assert(pUval != nullptr);
    if (!QuICC::Cuda::isDeviceMemory(pUval->data))
    {
        assert(!QuICC::Cuda::isDeviceMemory(pUmod->data));
        _ciface_quiccir_mul_const_buoyancy_f64_DCCSC3D_f64_DCCSC3D_cpu(obj, pUmod, pUval);
    }
    else
    {
        _ciface_quiccir_mul_const_buoyancy_f64_DCCSC3D_f64_DCCSC3D_gpu(obj, pUmod, pUval);
    }
    #else
    _ciface_quiccir_mul_const_buoyancy_f64_DCCSC3D_f64_DCCSC3D_cpu(obj, pUmod, pUval);
    #endif
}

/// @brief C Interface to MLIR for a mul const operator
/// @param obj
/// @param pRet
/// @param pU0
/// @param pU1
/// @param pU2
/// @param pV0
/// @param pV1
/// @param pV2
extern "C" void _ciface_quiccir_mul_const_transport_f64_DCCSC3D_f64_DCCSC3D(void* obj,
    ViewDescriptor<double, std::uint32_t, 3>* pRet,
    ViewDescriptor<double, std::uint32_t, 3>* pRhs)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_mul_const_transport_f64_DCCSC3D_f64_DCCSC3D\n";
    #endif
    assert(obj != nullptr);
    assert(pRet != nullptr);
    assert(pRhs != nullptr);
    // op
    using namespace QuICC::Slicewise::Cpu;
    using namespace QuICC::Slicewise;
    using T = QuICC::View::View<double, QuICC::View::DCCSC3D>;
    using op_t = Op<MulRFunctor<double>, T, T>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    ViewBase<std::uint32_t> pointers[rank];
    pointers[1] = ViewBase<std::uint32_t>(pRhs->pos, pRhs->posSize);
    ViewBase<std::uint32_t> indices[rank];
    indices[1] = ViewBase<std::uint32_t>(pRhs->coo, pRhs->cooSize);
    T viewRhs(pRhs->data, pRhs->dataSize, pRhs->dims, pointers, indices);
    T viewRet(pRet->data, pRet->dataSize, pRet->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewRet, viewRhs);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a binary mul const operator
/// gpu backend
/// @param obj
/// @param pRet
/// @param pRhs
extern "C" void _ciface_quiccir_mul_const_transport_f64_DCCSC3D_f64_DCCSC3D_gpu(void* obj,
    ViewDescriptor<double, std::uint32_t, 3>* pRet,
    ViewDescriptor<double, std::uint32_t, 3>* pRhs)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_mul_const_transport_f64_DCCSC3D_f64_DCCSC3D_gpu\n";
    #endif
    assert(obj != nullptr);
    assert(pRet != nullptr);
    assert(pRhs != nullptr);
    // op
    using namespace QuICC::Pointwise::Cuda;
    using namespace QuICC::Pointwise;
    using T = QuICC::View::View<double, QuICC::View::DCCSC3D>;
    using op_t = Op<MulRFunctor<double>, T, T>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    ViewBase<std::uint32_t> pointers[rank];
    pointers[1] = ViewBase<std::uint32_t>(pLhs->pos, pLhs->posSize);
    ViewBase<std::uint32_t> indices[rank];
    indices[1] = ViewBase<std::uint32_t>(pLhs->coo, pLhs->cooSize);
    T viewRhs(reinterpret_cast<cuda::double*>(pRhs->data), pRhs->dataSize, pRhs->dims, pointers, indices);
    T viewRet(reinterpret_cast<cuda::double*>(pRet->data), pRet->dataSize, pRet->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewRet, viewRhs);
};
#endif

/// @brief C Interface to MLIR for a mul const operator
/// @param op
/// @param uval
/// @param umod
extern "C" void _ciface_quiccir_mul_const_transport_f64_DCCSC3D_f64_DCCSC3D(void* obj, view3_cd_t* pUmod, view3_t* pUval)
{
    #ifdef QUICC_HAS_CUDA_BACKEND
    assert(pUval != nullptr);
    if (!QuICC::Cuda::isDeviceMemory(pUval->data))
    {
        assert(!QuICC::Cuda::isDeviceMemory(pUmod->data));
        _ciface_quiccir_mul_const_transport_f64_DCCSC3D_f64_DCCSC3D_cpu(obj, pUmod, pUval);
    }
    else
    {
        _ciface_quiccir_mul_const_transport_f64_DCCSC3D_f64_DCCSC3D_gpu(obj, pUmod, pUval);
    }
    #else
    _ciface_quiccir_mul_const_transport_f64_DCCSC3D_f64_DCCSC3D_cpu(obj, pUmod, pUval);
    #endif
}
