#include <iostream>
#include <complex>
#include <cassert>
#include <Quiccir-c/Utils.h>

#include "Graph/MlirShims.hpp"
#include "Graph/BackendsMap.hpp"
#include "View/View.hpp"

using namespace QuICC::Graph;


/// @brief C Interface to MLIR for a transpose operator
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_201_C_DCCSC3D_t_C_S1CLCSC3D_t(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_201_C_DCCSC3D_t_C_S1CLCSC3D_t\n";
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

extern "C" void _ciface_quiccir_transpose_201_C_DCCSC3D_t_C_DCCSC3D_t(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_201_C_DCCSC3D_t_C_DCCSC3D_t\n";
    #endif
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    assert(pIn->dataSize == pOut->dataSize);
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
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices);
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
};

extern "C" void _ciface_quiccir_alloc_transpose_201_C_DCCSC3D_t_C_S1CLCSC3D_t(view3_cd_t* pNewBuffer, const view3_cd_t* pProdBuffer)
{
    // This operation allocates for the serial transpose operator
    // therefore it assumes dense 3D tensors

    // Check slice dimensions according to permutation
    // note that the mangling comes from the MLIR
    // convention where the layers are leftmost
    assert(pNewBuffer->dims[2] == pProdBuffer->dims[0]);
    assert(pNewBuffer->dims[0] == pProdBuffer->dims[1]);
    assert(pNewBuffer->dims[1] == pProdBuffer->dims[2]);

    // Alloc meta for fully populated tensor
    // we will need to add a hashmap with counter
    // to know when we this was already allocated
    // and when we can deallocate the meta data
    auto meta = denseTransposePtrAndIdx<C_DCCSC3D_t, C_S1CLCSC3D_t>(
        {pNewBuffer->dims[0], pNewBuffer->dims[1], pNewBuffer->dims[2]});

    pNewBuffer->posSize = meta.ptr.size();
    std::size_t sizeByte = sizeof(std::uint32_t) * pNewBuffer->posSize;
    pNewBuffer->pos = reinterpret_cast<std::uint32_t*>(::operator new(sizeByte, static_cast<std::align_val_t>(sizeof(std::uint32_t))));
    pNewBuffer->cooSize = meta.idx.size();
    sizeByte = sizeof(std::uint32_t) * pNewBuffer->cooSize;
    pNewBuffer->coo = reinterpret_cast<std::uint32_t*>(::operator new(sizeByte, static_cast<std::align_val_t>(sizeof(std::uint32_t))));

    // Populate meta for fully populated tensor
    for (std::size_t i = 1; i < pNewBuffer->posSize; ++i) {
        pNewBuffer->pos[i] = meta.ptr[i];
    }
    for (std::size_t i = 0; i < pNewBuffer->cooSize; ++i) {
        pNewBuffer->coo[i] = meta.idx[i];
    }

    // Alloc buffer
    pNewBuffer->dataSize = pNewBuffer->dims[0] * pNewBuffer->cooSize;
    sizeByte = sizeof(std::complex<double>) * pNewBuffer->dataSize;
    pNewBuffer->data = reinterpret_cast<std::complex<double>*>(::operator new(sizeByte, static_cast<std::align_val_t>(sizeof(std::complex<double>))));
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_transpose_201_C_DCCSC3D_t_C_S1CLCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};


extern "C" void _ciface_quiccir_alloc_transpose_201_C_DCCSC3D_t_C_DCCSC3D_t(view3_cd_t* pNewBuffer, const view3_cd_t* pProdBuffer)
{
    // This operation allocates for the serial transpose operator
    // therefore it assumes dense 3D tensors

    // Check slice dimensions according to permutation
    // note that the mangling comes from the MLIR
    // convention where the layers are leftmost
    assert(pNewBuffer->dims[2] == pProdBuffer->dims[0]);
    assert(pNewBuffer->dims[0] == pProdBuffer->dims[1]);
    assert(pNewBuffer->dims[1] == pProdBuffer->dims[2]);

    // Alloc meta for fully populated tensor
    // we will need to add a hashmap with counter
    // to know when we can deallocate the meta data
    pNewBuffer->posSize = pNewBuffer->dims[2]+1;
    std::size_t sizeByte = sizeof(std::uint32_t) * pNewBuffer->posSize;
    pNewBuffer->pos = reinterpret_cast<std::uint32_t*>(::operator new(sizeByte, static_cast<std::align_val_t>(sizeof(std::uint32_t))));

    pNewBuffer->cooSize = pNewBuffer->dims[1]*pNewBuffer->dims[2];
    sizeByte = sizeof(std::uint32_t) * pNewBuffer->posSize;
    pNewBuffer->coo = reinterpret_cast<std::uint32_t*>(::operator new(sizeByte, static_cast<std::align_val_t>(sizeof(std::uint32_t))));
    // Populate meta for fully populated tensor
    pNewBuffer->pos[0] = 0;
    for (std::size_t i = 1; i < pNewBuffer->posSize; ++i) {
        pNewBuffer->pos[i] = pNewBuffer->pos[i-1]+pNewBuffer->dims[1];
    }
    for (std::size_t i = 0; i < pNewBuffer->cooSize; ++i) {
        pNewBuffer->coo[i] = i % pNewBuffer->dims[1];
    }

    // Alloc buffer
    pNewBuffer->dataSize = pNewBuffer->dims[0] * pNewBuffer->cooSize;
    sizeByte = sizeof(std::complex<double>) * pNewBuffer->dataSize;
    pNewBuffer->data = reinterpret_cast<std::complex<double>*>(::operator new(sizeByte, static_cast<std::align_val_t>(sizeof(std::complex<double>))));
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_transpose_201_C_DCCSC3D_t_C_DCCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};
