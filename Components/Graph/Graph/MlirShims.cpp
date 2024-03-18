#include <iostream>
#include <Quiccir-c/Utils.h>
#include <complex>

#include "Graph/MlirShims.hpp"
#include "Graph/BackendsMap.hpp"
#include "View/View.hpp"

using namespace QuICC::Graph;

/// @brief C Interface to MLIR for a fr prj operator
/// @param op
/// @param uval
/// @param umod
extern "C" void _ciface_quiccir_fr_prj_R_DCCSC3D_t_C_DCCSC3D_t(void* obj, view3_t* pUval, view3_cd_t* pUmod)
{
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

/// @brief C Interface to MLIR for a fr int operator
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_fr_int_C_DCCSC3D_t_R_DCCSC3D_t(void* obj, view3_cd_t* pUmod, view3_t* pUval)
{
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


/// @brief C Interface to MLIR for an allocator
/// The name fully describes how the buffer needs to be allocated
/// since it includes the producer op and buffer.
/// @param pNewBuffer buffer to be allocated
/// @param pProdBuffer producer buffer
extern "C" void _ciface_quiccir_alloc_fr_prj_R_DCCSC3D_t_C_DCCSC3D_t(view3_t* pNewBuffer, view3_cd_t* pProdBuffer)
{
    // Slice (phys = mods*2-1)
    assert(pNewBuffer->dims[0] == pProdBuffer->dims[0]*2-1);
    assert(pNewBuffer->dims[1] == pProdBuffer->dims[1]);
    // Layers
    assert(pNewBuffer->dims[2] == pProdBuffer->dims[2]);
    // reuse meta
    pNewBuffer->pos = pProdBuffer->pos;
    pNewBuffer->posSize = pProdBuffer->posSize;
    pNewBuffer->coo = pProdBuffer->coo;
    pNewBuffer->cooSize = pProdBuffer->cooSize;
    // alloc buffer
    pNewBuffer->dataSize = pNewBuffer->dims[0] * pNewBuffer->cooSize;
    std::size_t sizeByte = sizeof(double) * pNewBuffer->dataSize;
    pNewBuffer->data = reinterpret_cast<double*>(::operator new(sizeByte, static_cast<std::align_val_t>(sizeof(double))));
};

/// @brief C Interface to MLIR for a  deallocator
/// @param umod
extern "C" void _ciface_quiccir_dealloc_R_DCCSC3D_t(view3_t* pBuffer)
{
    // meta
    pBuffer->coo = nullptr;
    pBuffer->cooSize = 0;
    pBuffer->pos = nullptr;
    pBuffer->posSize = 0;
    // dealloc
    std::size_t sizeByte = sizeof(double) * pBuffer->dataSize;
    ::operator delete(pBuffer->data, sizeByte, static_cast<std::align_val_t>(sizeof(double)));
    pBuffer->dataSize = 0;
};