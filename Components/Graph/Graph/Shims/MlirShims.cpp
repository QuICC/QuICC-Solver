#include <iostream>
#include <complex>
#include <cassert>

#include "Graph/Shims/MlirShims.hpp"
#include "Graph/BackendsMap.hpp"
#include "Graph/Types.hpp"

using namespace QuICC::Graph;

/// @brief C Interface to MLIR for a binary add operator
/// @param obj
/// @param pRet
/// @param pLhs
/// @param pRhs
extern "C" void _ciface_quiccir_add_C_DCCSC3D_t_C_DCCSC3D_t_C_DCCSC3D_t(void* obj,
    view3_cd_t* pRet, view3_cd_t* pLhs, view3_cd_t* pRhs)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_add_C_DCCSC3D_t_C_DCCSC3D_t_C_DCCSC3D_t\n";
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

/// @brief C Interface to MLIR for a binary sub operator
/// @param obj
/// @param pRet
/// @param pLhs
/// @param pRhs
extern "C" void _ciface_quiccir_sub_C_DCCSC3D_t_C_DCCSC3D_t_C_DCCSC3D_t(void* obj,
    view3_cd_t* pRet, view3_cd_t* pLhs, view3_cd_t* pRhs)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_sub_C_DCCSC3D_t_C_DCCSC3D_t_C_DCCSC3D_t\n";
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

/// @brief C Interface to MLIR for an allocator
/// @param pNewBuffer
/// @param pProdBuffer
extern "C" void _ciface_quiccir_alloc_add_C_DCCSC3D_t_C_DCCSC3D_t(view3_cd_t* pNewBuffer, view3_cd_t* pProdBuffer)
{
    // NewBuffer is identical to producer
    assert(pNewBuffer->dims[0] == pProdBuffer->dims[0]);
    assert(pNewBuffer->dims[1] == pProdBuffer->dims[1]);
    assert(pNewBuffer->dims[2] == pProdBuffer->dims[2]);
    // Reuse meta
    assert(pProdBuffer->pos != nullptr);
    assert(pProdBuffer->coo != nullptr);
    pNewBuffer->pos = pProdBuffer->pos;
    pNewBuffer->posSize = pProdBuffer->posSize;
    pNewBuffer->coo = pProdBuffer->coo;
    pNewBuffer->cooSize = pProdBuffer->cooSize;
    // Alloc buffer
    pNewBuffer->dataSize = pProdBuffer->dataSize;
    std::size_t sizeByte = sizeof(std::complex<double>) * pNewBuffer->dataSize;
    pNewBuffer->data = reinterpret_cast<std::complex<double>*>(::operator new(sizeByte, static_cast<std::align_val_t>(sizeof(std::complex<double>))));
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_add_C_DCCSC3D_t_C_DCCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};

/// @brief C Interface to MLIR for an allocator
/// @param pNewBuffer
/// @param pProdBuffer
extern "C" void _ciface_quiccir_alloc_sub_C_DCCSC3D_t_C_DCCSC3D_t(view3_cd_t* pNewBuffer, view3_cd_t* pProdBuffer)
{
    // NewBuffer is identical to producer
    assert(pNewBuffer->dims[0] == pProdBuffer->dims[0]);
    assert(pNewBuffer->dims[1] == pProdBuffer->dims[1]);
    assert(pNewBuffer->dims[2] == pProdBuffer->dims[2]);
    // Reuse meta
    assert(pProdBuffer->pos != nullptr);
    assert(pProdBuffer->coo != nullptr);
    pNewBuffer->pos = pProdBuffer->pos;
    pNewBuffer->posSize = pProdBuffer->posSize;
    pNewBuffer->coo = pProdBuffer->coo;
    pNewBuffer->cooSize = pProdBuffer->cooSize;
    // Alloc buffer
    pNewBuffer->dataSize = pProdBuffer->dataSize;
    std::size_t sizeByte = sizeof(std::complex<double>) * pNewBuffer->dataSize;
    pNewBuffer->data = reinterpret_cast<std::complex<double>*>(::operator new(sizeByte, static_cast<std::align_val_t>(sizeof(std::complex<double>))));
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_sub_C_DCCSC3D_t_C_DCCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};
