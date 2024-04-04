#include <iostream>
#include <complex>
#include <cassert>
#include <Quiccir-c/Utils.h>

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

/// @brief C Interface to MLIR for a fr int operator
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_fr_int_C_DCCSC3D_t_R_DCCSC3D_t(void* obj, view3_cd_t* pUmod, view3_t* pUval)
{
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

/// @brief C Interface to MLIR for a al int operator
/// @param op
/// @param uval
/// @param umod
extern "C" void _ciface_quiccir_al_int_C_S1CLCSC3D_t_C_DCCSC3D_t(void* obj, view3_cd_t* pUval, view3_cd_t* pUmod)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_al_int_C_S1CLCSC3D_t_C_DCCSC3D_t\n";
    #endif
    assert(obj != nullptr);
    assert(pUval != nullptr);
    assert(pUmod != nullptr);
    // op
    using namespace QuICC::Transform::Quadrature;
    using Tin = C_DCCSC3D_t;
    using Tout = C_S1CLCSC3D_t;
    using Top = QuICC::View::View<double, QuICC::View::S1CLCSC3DJIK>;
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
    if (cl->getOp().data() == nullptr) {
        // set up operator
        constexpr size_t rank = 3;
        /// dim 0 - L  - harmonic degree
        /// dim 1 - Nl - longitudinal points
        /// dim 2 - M  - harmonic order
        std::array<std::uint32_t, rank> dims {pUval->dims[0], pUmod->dims[0], pUval->dims[2]};
        std::vector<std::uint32_t> layers;
        // check for populated layers
        for (std::size_t i = 0; i < pUmod->posSize-1; ++i) {
            if (pUmod->pos[i+1] - pUmod->pos[i] > 0) {
                layers.push_back(i);
            }
        }
        cl->allocOp(dims, layers);
    }
    // call
    cl->apply(viewVal, viewMod);
};

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
    assert(cl->getOp().data() == nullptr);
    // call
    cl->apply(viewVal, viewMod);
};

/// @brief C Interface to MLIR for a binary add operator
/// @param op
/// @param uval
/// @param umod
extern "C" void _ciface_quiccir_add_C_DCCSC3D_t_C_DCCSC3D_t_C_DCCSC3D_t(void* obj,
    view3_cd_t* pRet, view3_cd_t* pLhs, view3_cd_t* pRhs)
{
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
/// @param op
/// @param uval
/// @param umod
extern "C" void _ciface_quiccir_sub_C_DCCSC3D_t_C_DCCSC3D_t_C_DCCSC3D_t(void* obj,
    view3_cd_t* pRet, view3_cd_t* pLhs, view3_cd_t* pRhs)
{
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



/// @brief C Interface to MLIR for a transpose operator
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_120_C_DCCSC3D_t_C_S1CLCSC3D_t(void* obj, view3_cd_t*, view3_cd_t*)
{
   std::cout << "missing transpose op shim\n";
};

extern "C" void _ciface_quiccir_transpose_021_C_DCCSC3D_t_C_DCCSC3D_t(void* obj, view3_cd_t* pOut, view3_cd_t* pIn)
{
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    // op
    using namespace QuICC::Transpose::Cpu;
    using namespace QuICC::Transpose;
    using Tin = C_DCCSC3D_t;
    using Tout = C_DCCSC3D_t;
    using op_t = Op<Tout, Tin, p021_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    // not used, not setting up
    ViewBase<std::uint32_t> pointers[rank];
    // not used, not setting up
    ViewBase<std::uint32_t> indices[rank];
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices);
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_021_C_DCCSC3D_t_C_DCCSC3D_t\n";
    #endif
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
    // Reuse meta
    assert(pProdBuffer->pos != nullptr);
    assert(pProdBuffer->coo != nullptr);
    pNewBuffer->pos = pProdBuffer->pos;
    pNewBuffer->posSize = pProdBuffer->posSize;
    pNewBuffer->coo = pProdBuffer->coo;
    pNewBuffer->cooSize = pProdBuffer->cooSize;
    // Alloc buffer
    pNewBuffer->dataSize = pNewBuffer->dims[0] * pNewBuffer->cooSize;
    std::size_t sizeByte = sizeof(double) * pNewBuffer->dataSize;
    pNewBuffer->data = reinterpret_cast<double*>(::operator new(sizeByte, static_cast<std::align_val_t>(sizeof(double))));
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_fr_prj_R_DCCSC3D_t_C_DCCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};

extern "C" void _ciface_quiccir_alloc_fr_int_C_DCCSC3D_t_R_DCCSC3D_t(view3_cd_t* pNewBuffer, view3_t* pProdBuffer)
{
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
    std::size_t sizeByte = sizeof(std::complex<double>) * pNewBuffer->dataSize;
    pNewBuffer->data = reinterpret_cast<std::complex<double>*>(::operator new(sizeByte, static_cast<std::align_val_t>(sizeof(std::complex<double>))));
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_fr_int_C_DCCSC3D_t_R_DCCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};

extern "C" void _ciface_quiccir_alloc_al_int_C_S1CLCSC3D_t_C_DCCSC3D_t(view3_cd_t* pNewBuffer, view3_cd_t* pProdBuffer)
{
    // Check slice dimensions
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
    std::size_t cumSliceSize = 0;
    for (std::size_t i = 0; i < pNewBuffer->posSize - 1; ++i) {
        auto width = pNewBuffer->pos[i+1] - pNewBuffer->pos[i];
        auto height = pNewBuffer->dims[0] - i;
        cumSliceSize += height * width;
    }
    pNewBuffer->dataSize = cumSliceSize;
    std::size_t sizeByte = sizeof(std::complex<double>) * pNewBuffer->dataSize;
    pNewBuffer->data = reinterpret_cast<std::complex<double>*>(::operator new(sizeByte, static_cast<std::align_val_t>(sizeof(std::complex<double>))));
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_al_int_C_S1CLCSC3D_t_C_DCCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};

extern "C" void _ciface_quiccir_alloc_jw_int_C_DCCSC3D_t_C_DCCSC3D_t(view3_cd_t* pNewBuffer, view3_cd_t* pProdBuffer)
{
    std::cout << "missing alloc operator shim: _ciface_quiccir_alloc_jw_int_C_DCCSC3D_t_C_DCCSC3D_t\n";
    pNewBuffer->pos = nullptr;
    pNewBuffer->posSize = 0;
    pNewBuffer->coo = nullptr;
    pNewBuffer->cooSize = 0;
};

extern "C" void _ciface_quiccir_alloc_sub_C_DCCSC3D_t_C_DCCSC3D_t(view3_cd_t* pNewBuffer, view3_cd_t* pProdBuffer)
{
    std::cout << "missing alloc operator shim: _ciface_quiccir_alloc_sub_C_DCCSC3D_t_C_DCCSC3D_t\n";
    pNewBuffer->pos = nullptr;
    pNewBuffer->posSize = 0;
    pNewBuffer->coo = nullptr;
    pNewBuffer->cooSize = 0;
};

extern "C" void _ciface_quiccir_alloc_transpose_021_C_DCCSC3D_t_C_DCCSC3D_t(view3_cd_t* pNewBuffer, view3_cd_t* pProdBuffer)
{
    // This operation allocates for the serial transpose operator
    // therefore it assumes dense 3D tensors

    // Check slice dimensions according to permutation
    // note that the mangling comes from the MLIR
    // convention where the layers are leftmost
    assert(pNewBuffer->dims[1] == pProdBuffer->dims[0]);
    assert(pNewBuffer->dims[0] == pProdBuffer->dims[1]);
    assert(pNewBuffer->dims[2] == pProdBuffer->dims[2]);

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
    std::cout << "_ciface_quiccir_alloc_transpose_021_C_DCCSC3D_t_C_DCCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};

extern "C" void _ciface_quiccir_alloc_transpose_120_C_DCCSC3D_t_C_S1CLCSC3D_t(view3_cd_t* pNewBuffer, view3_cd_t* pProdBuffer)
{
    std::cout << "missing alloc operator shim: _ciface_quiccir_alloc_transpose_120_C_DCCSC3D_t_C_S1CLCSC3D_t\n";
    pNewBuffer->pos = nullptr;
    pNewBuffer->posSize = 0;
    pNewBuffer->coo = nullptr;
    pNewBuffer->cooSize = 0;
};


/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_R_DCCSC3D_t(view3_t* pBuffer)
{
    // meta
    pBuffer->coo = nullptr;
    pBuffer->cooSize = 0;
    pBuffer->pos = nullptr;
    pBuffer->posSize = 0;
    // dealloc
    assert(pBuffer->data != nullptr);
    std::size_t sizeByte = sizeof(double) * pBuffer->dataSize;
    ::operator delete(pBuffer->data, sizeByte, static_cast<std::align_val_t>(sizeof(double)));
    pBuffer->dataSize = 0;
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_R_DCCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_C_DCCSC3D_t(view3_t* pBuffer)
{
    // meta
    pBuffer->coo = nullptr;
    pBuffer->cooSize = 0;
    pBuffer->pos = nullptr;
    pBuffer->posSize = 0;
    // dealloc
    assert(pBuffer->data != nullptr);
    std::size_t sizeByte = sizeof(std::complex<double>) * pBuffer->dataSize;
    ::operator delete(pBuffer->data, sizeByte, static_cast<std::align_val_t>(sizeof(std::complex<double>)));
    pBuffer->dataSize = 0;
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_C_DCCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};

/// @brief C Interface to MLIR for a deallocator
/// @param pBuffer
extern "C" void _ciface_quiccir_dealloc_C_S1CLCSC3D_t(view3_t* pBuffer)
{
    // meta
    pBuffer->coo = nullptr;
    pBuffer->cooSize = 0;
    pBuffer->pos = nullptr;
    pBuffer->posSize = 0;
    // dealloc
    assert(pBuffer->data != nullptr);
    std::size_t sizeByte = sizeof(std::complex<double>) * pBuffer->dataSize;
    ::operator delete(pBuffer->data, sizeByte, static_cast<std::align_val_t>(sizeof(std::complex<double>)));
    pBuffer->dataSize = 0;
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_dealloc_C_S1CLCSC3D_t, bytes: " << sizeByte << '\n';
    #endif
};
