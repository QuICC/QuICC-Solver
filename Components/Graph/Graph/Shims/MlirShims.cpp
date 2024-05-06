#include <iostream>
#include <complex>
#include <cassert>

#include "Graph/Shims/MlirShims.hpp"
#include "Graph/BackendsMap.hpp"
#include "Graph/Types.hpp"


namespace QuICC {
namespace Graph {

/// @brief compute pointers and index meta data for fully populated
/// tensor as ouput of a Transpose operation between AL and JW
/// @tparam C_DCCSC3D_t new buffer type
/// @tparam C_S1CLCSC3D_t producer type
/// @param dims output dimensions
/// @return
template<>
ptrAndIdx denseTransposePtrAndIdx<C_DCCSC3D_t, C_S1CLCSC3D_t>(
    const std::array<std::uint32_t, 3> dims)
{
    ptrAndIdx ret;
    std::uint32_t K = dims[1];
    std::uint32_t I = dims[2];

    // row width (with jki) K - ...
    std::vector<std::uint32_t> kLess(I);
    kLess[I-1] = K-1;
    for (std::size_t i = I-1; i > 0; --i)
    {
        if (kLess[i] > 0)
        {
            kLess[i-1] = kLess[i] - 1;
        }
        else
        {
            kLess[i-1] = 0;
        }
    }

    std::size_t layWidthCum = 0;
    ret.ptr.resize(I+1);
    ret.ptr[0] = 0;
    for (std::size_t i = 1; i < I+1; ++i)
    {
        std::uint32_t width = K-kLess[i-1];
        ret.ptr[i] = ret.ptr[i-1] + width;
        layWidthCum += width;
    }

    std::size_t layerIdx = 0;
    ret.idx.resize(layWidthCum);
    for (std::size_t l = 0; l < ret.ptr.size() - 1; ++l)
    {
        auto layerSize = ret.ptr[l+1] - ret.ptr[l];
        for (std::size_t i = 0; i < layerSize; ++i)
        {
            ret.idx[layerIdx+i] = i;
        }
        layerIdx += layerSize;
    }
    return ret;
}

template<>
ptrAndIdx denseTransposePtrAndIdx<C_S1CLCSC3D_t, C_DCCSC3D_t>(const std::array<std::uint32_t, 3> dims)
{
    ptrAndIdx ret;
    std::uint32_t J = dims[1];
    std::uint32_t K = dims[2];

    // fully populated, each layer has all columns
    ret.ptr.resize(K+1);
    ret.ptr[0] = 0;
    for (std::size_t k = 1; k < ret.ptr.size(); ++k)
    {
        ret.ptr[k] = ret.ptr[k-1] + J;
    }

    ret.idx.resize(J*K);
    for (std::size_t i = 1; i < ret.idx.size(); ++i)
    {
        ret.idx[i] = i % J;
    }
    return ret;
}

} // namespace
} // namespace


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
