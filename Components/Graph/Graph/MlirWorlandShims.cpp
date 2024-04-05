#include <iostream>
#include <complex>
#include <cassert>
#include <Quiccir-c/Utils.h>

#include "Graph/MlirShims.hpp"
#include "Graph/BackendsMap.hpp"
#include "View/View.hpp"

using namespace QuICC::Graph;

/// @brief C Interface to MLIR for a jw int operator
/// @param op
/// @param uval
/// @param umod
extern "C" void _ciface_quiccir_jw_int_C_DCCSC3D_t_C_DCCSC3D_t(void* obj, view3_cd_t* pUval, view3_cd_t* pUmod)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_jw_int_C_DCCSC3D_t_C_DCCSC3D_t\n";
    #endif
    assert(obj != nullptr);
    assert(pUval != nullptr);
    assert(pUmod != nullptr);
    // op
    using namespace QuICC::Transform::Quadrature;
    using Tin = C_DCCSC3D_t;
    using Tout = C_DCCSC3D_t;
    using Top = QuICC::View::View<double, QuICC::View::CSL3DJIK>;
    using backend_t = Cpu::ImplOp<Tout, Tin, Top>;
    using op_t = Op<Tout, Tin, Top, backend_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    ViewBase<std::uint32_t> pointers[rank];
    pointers[1] = ViewBase<std::uint32_t>(pUmod->pos, pUmod->posSize);
    ViewBase<std::uint32_t> indices[rank];
    indices[1] = ViewBase<std::uint32_t>(pUmod->coo, pUmod->cooSize);
    Tin viewMod(pUmod->data, pUmod->dataSize, pUmod->dims, pointers, indices);
    Tout viewVal(pUval->data, pUval->dataSize, pUval->dims, pointers, indices);
    // Check that op was set up
    auto cl = reinterpret_cast<op_t*>(obj);
    assert(cl->getOp().data() != nullptr);
    // call
    cl->apply(viewVal, viewMod);
};

extern "C" void _ciface_quiccir_alloc_jw_int_C_DCCSC3D_t_C_DCCSC3D_t(view3_cd_t* pNewBuffer, view3_cd_t* pProdBuffer)
{
    // Integrator, NewBuffer is in modal space
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
    std::size_t sizeByte = sizeof(std::complex<double>) * pNewBuffer->dataSize;
    pNewBuffer->data = reinterpret_cast<std::complex<double>*>(::operator new(sizeByte, static_cast<std::align_val_t>(sizeof(std::complex<double>))));
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_jw_int_C_DCCSC3D_t_C_DCCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};
