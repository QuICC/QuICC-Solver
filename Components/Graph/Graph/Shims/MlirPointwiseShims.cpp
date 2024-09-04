#include <iostream>
#include <complex>
#include <cassert>

#include "Graph/Shims/MlirShims.hpp"
#include "ViewOps/Pointwise/Pointwise.hpp"
#include "ViewOps/Pointwise/Functors.hpp"
#include "Graph/Types.hpp"


using namespace QuICC::Graph;

/// @brief C Interface to MLIR for a binary add operator
/// @param obj
/// @param pRet
/// @param pLhs
/// @param pRhs
extern "C" void _ciface_quiccir_add_complexf64_DCCSC3D_complexf64_DCCSC3D_complexf64_DCCSC3D(void* obj,
    view3_cd_t* pRet, view3_cd_t* pLhs, view3_cd_t* pRhs)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_add_complexf64_DCCSC3D_complexf64_DCCSC3D_complexf64_DCCSC3D\n";
    #endif
    assert(obj != nullptr);
    assert(pRet != nullptr);
    assert(pLhs != nullptr);
    assert(pRhs != nullptr);
    // op
    using namespace QuICC::Pointwise::Cpu;
    using namespace QuICC::Pointwise;
    using T = C_DCCSC3D_t;
    using op_t = Op<AddFunctor<std::complex<double>>, T, T, T>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    ViewBase<std::uint32_t> pointers[rank];
    pointers[1] = ViewBase<std::uint32_t>(pLhs->pos, pLhs->posSize);
    ViewBase<std::uint32_t> indices[rank];
    indices[1] = ViewBase<std::uint32_t>(pLhs->coo, pLhs->cooSize);
    T viewLhs(pLhs->data, pLhs->dataSize, pLhs->dims, pointers, indices);
    T viewRhs(pRhs->data, pRhs->dataSize, pRhs->dims, pointers, indices);
    T viewRet(pRet->data, pRet->dataSize, pRet->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewRet, viewLhs, viewRhs);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a binary add operator
/// gpu backend
/// @param obj
/// @param pRet
/// @param pLhs
/// @param pRhs
extern "C" void _ciface_quiccir_add_complexf64_DCCSC3DJIK_complexf64_DCCSC3DJIK_complexf64_DCCSC3DJIK(void* obj,
    view3_cd_t* pRet, view3_cd_t* pLhs, view3_cd_t* pRhs)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_add_complexf64_DCCSC3DJIK_complexf64_DCCSC3DJIK_complexf64_DCCSC3DJIK\n";
    #endif
    assert(obj != nullptr);
    assert(pRet != nullptr);
    assert(pLhs != nullptr);
    assert(pRhs != nullptr);
    // op
    using namespace QuICC::Pointwise::Cuda;
    using namespace QuICC::Pointwise;
    using T = Ccuda_DCCSC3DJIK_t;
    using op_t = Op<AddFunctor<typename T::ScalarType>, T, T, T>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    ViewBase<std::uint32_t> pointers[rank];
    pointers[1] = ViewBase<std::uint32_t>(pLhs->pos, pLhs->posSize);
    ViewBase<std::uint32_t> indices[rank];
    indices[1] = ViewBase<std::uint32_t>(pLhs->coo, pLhs->cooSize);
    T viewLhs(reinterpret_cast<cuda::std::complex<double>*>(pLhs->data), pLhs->dataSize, pLhs->dims, pointers, indices);
    T viewRhs(reinterpret_cast<cuda::std::complex<double>*>(pRhs->data), pRhs->dataSize, pRhs->dims, pointers, indices);
    T viewRet(reinterpret_cast<cuda::std::complex<double>*>(pRet->data), pRet->dataSize, pRet->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewRet, viewLhs, viewRhs);
};
#endif

/// @brief C Interface to MLIR for a binary sub operator
/// @param obj
/// @param pRet
/// @param pLhs
/// @param pRhs
extern "C" void _ciface_quiccir_sub_complexf64_DCCSC3D_complexf64_DCCSC3D_complexf64_DCCSC3D(void* obj,
    view3_cd_t* pRet, view3_cd_t* pLhs, view3_cd_t* pRhs)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_sub_complexf64_DCCSC3D_complexf64_DCCSC3D_complexf64_DCCSC3D\n";
    #endif
    assert(obj != nullptr);
    assert(pRet != nullptr);
    assert(pLhs != nullptr);
    assert(pRhs != nullptr);
    // op
    using namespace QuICC::Pointwise::Cpu;
    using namespace QuICC::Pointwise;
    using T = C_DCCSC3D_t;
    using op_t = Op<SubFunctor<std::complex<double>>, T, T, T>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    ViewBase<std::uint32_t> pointers[rank];
    pointers[1] = ViewBase<std::uint32_t>(pLhs->pos, pLhs->posSize);
    ViewBase<std::uint32_t> indices[rank];
    indices[1] = ViewBase<std::uint32_t>(pLhs->coo, pLhs->cooSize);
    T viewLhs(pLhs->data, pLhs->dataSize, pLhs->dims, pointers, indices);
    T viewRhs(pRhs->data, pRhs->dataSize, pRhs->dims, pointers, indices);
    T viewRet(pRet->data, pRet->dataSize, pRet->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewRet, viewLhs, viewRhs);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a binary sub operator
/// gpu backend
/// @param obj
/// @param pRet
/// @param pLhs
/// @param pRhs
extern "C" void _ciface_quiccir_sub_complexf64_DCCSC3DJIK_complexf64_DCCSC3DJIK_complexf64_DCCSC3DJIK(void* obj,
    view3_cd_t* pRet, view3_cd_t* pLhs, view3_cd_t* pRhs)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_sub_complexf64_DCCSC3DJIK_complexf64_DCCSC3DJIK_complexf64_DCCSC3DJIK\n";
    #endif
    assert(obj != nullptr);
    assert(pRet != nullptr);
    assert(pLhs != nullptr);
    assert(pRhs != nullptr);
    // op
    using namespace QuICC::Pointwise::Cuda;
    using namespace QuICC::Pointwise;
    using T = Ccuda_DCCSC3DJIK_t;
    using op_t = Op<SubFunctor<typename T::ScalarType>, T, T, T>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    ViewBase<std::uint32_t> pointers[rank];
    pointers[1] = ViewBase<std::uint32_t>(pLhs->pos, pLhs->posSize);
    ViewBase<std::uint32_t> indices[rank];
    indices[1] = ViewBase<std::uint32_t>(pLhs->coo, pLhs->cooSize);
    T viewLhs(reinterpret_cast<cuda::std::complex<double>*>(pLhs->data), pLhs->dataSize, pLhs->dims, pointers, indices);
    T viewRhs(reinterpret_cast<cuda::std::complex<double>*>(pRhs->data), pRhs->dataSize, pRhs->dims, pointers, indices);
    T viewRet(reinterpret_cast<cuda::std::complex<double>*>(pRet->data), pRet->dataSize, pRet->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewRet, viewLhs, viewRhs);
};
#endif
