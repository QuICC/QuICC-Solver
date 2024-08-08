#include <iostream>
#include <complex>
#include <cassert>

#include "Graph/Shims/MlirShims.hpp"
#include "Graph/BackendsMap.hpp"
#include "Graph/Types.hpp"

#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandRule.hpp"
#include "ViewOps/Worland/Tags.hpp"
#include "ViewOps/Worland/Builder.hpp"
#include "DenseSM/Worland/Operator.hpp"

using namespace QuICC::Graph;

/// @brief C Interface to MLIR for a jw int operator
/// column major data, cpu operators
/// @param op
/// @param uval
/// @param umod
extern "C" void _ciface_quiccir_jw_int_complexf64_DCCSC3D_complexf64_DCCSC3D(void* obj, view3_cd_t* pUmod, view3_cd_t* pUval)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_jw_int_complexf64_DCCSC3D_complexf64_DCCSC3D\n";
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
    /// \todo stronger check?
    assert(pUmod->posSize == pUval->posSize);
    assert(pUmod->cooSize == pUval->cooSize);
    Tin viewVal(pUval->data, pUval->dataSize, pUval->dims, pointers, indices);
    Tout viewMod(pUmod->data, pUmod->dataSize, pUmod->dims, pointers, indices);
    // Check that op was set up
    auto cl = reinterpret_cast<op_t*>(obj);
    if (cl->getOp().data() == nullptr)
    {
        // dim 0 - N  - radial modes
        // dim 1 - Nr - radial points
        // dim 2 - L  - harmonic degree
        std::array<std::uint32_t, rank> dims {pUmod->dims[0], pUval->dims[0], pUmod->dims[2]};
        std::vector<std::uint32_t> layers;
        for (std::size_t i = 0; i < dims[2]; ++i) {
            layers.push_back(i);
        }
        cl->allocOp(dims, layers);
        /// Set grid \todo set once per operator kind
        ::QuICC::Internal::Array igrid;
        ::QuICC::Internal::Array iweights;
        ::QuICC::Polynomial::Quadrature::WorlandRule quad;
        quad.computeQuadrature(igrid, iweights, pUval->dims[0]);
        // Populate op
        auto opView = cl->getOp();
        using namespace QuICC::Polynomial::Worland;
        using namespace QuICC::Transform::Worland;
        Builder<Top, QuICC::DenseSM::Worland::Operator<Wnl>, fwd_t> tBuilder;
        tBuilder.compute(opView, igrid, iweights);
    }
    assert(cl->getOp().data() != nullptr);
    // call
    cl->apply(viewMod, viewVal);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a jw int operator
/// row major data, gpu operators
/// @param op
/// @param uval
/// @param umod
extern "C" void _ciface_quiccir_jw_int_complexf64_DCCSC3DJIK_complexf64_DCCSC3DJIK(void* obj, view3_cd_t* pUmod, view3_cd_t* pUval)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_jw_int_complexf64_DCCSC3DJIK_complexf64_DCCSC3DJIK\n";
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
    using namespace QuICC::Transform::Quadrature;
    using Tin = C_DCCSC3DJIK_t;
    using Tout = C_DCCSC3DJIK_t;
    using Top = QuICC::View::View<double, QuICC::View::CSL3D>;
    using backend_t = Cuda::ImplOp<Tout, Tin, Top>;
    using op_t = Op<Tout, Tin, Top, backend_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    ViewBase<std::uint32_t> pointers[rank];
    pointers[1] = ViewBase<std::uint32_t>(pUmod->pos, pUmod->posSize);
    ViewBase<std::uint32_t> indices[rank];
    indices[1] = ViewBase<std::uint32_t>(pUmod->coo, pUmod->cooSize);
    /// \todo stronger check?
    assert(pUmod->posSize == pUval->posSize);
    assert(pUmod->cooSize == pUval->cooSize);
    Tin viewVal(pUval->data, pUval->dataSize, pUval->dims, pointers, indices);
    Tout viewMod(pUmod->data, pUmod->dataSize, pUmod->dims, pointers, indices);
    // Check that op was set up
    auto cl = reinterpret_cast<op_t*>(obj);
    if (cl->getOp().data() == nullptr)
    {
        // dim 0 - N  - radial modes
        // dim 1 - Nr - radial points
        // dim 2 - L  - harmonic degree
        std::array<std::uint32_t, rank> dims {pUmod->dims[0], pUval->dims[0], pUmod->dims[2]};
        std::vector<std::uint32_t> layers;
        for (std::size_t i = 0; i < dims[2]; ++i) {
            layers.push_back(i);
        }
        cl->allocOp(dims, layers);
        /// Set grid \todo set once per operator kind
        ::QuICC::Internal::Array igrid;
        ::QuICC::Internal::Array iweights;
        ::QuICC::Polynomial::Quadrature::WorlandRule quad;
        quad.computeQuadrature(igrid, iweights, pUval->dims[0]);
        // Populate op
        auto opView = cl->getOp();
        using namespace QuICC::Polynomial::Worland;
        using namespace QuICC::Transform::Worland;
        Builder<Top, QuICC::DenseSM::Worland::Operator<Wnl>, fwd_t> tBuilder;
        tBuilder.compute(opView, igrid, iweights);
    }
    assert(cl->getOp().data() != nullptr);
    // call
    cl->apply(viewMod, viewVal);
};
#endif

/// @brief C Interface to MLIR for a jw prj operator
/// column major data, cpu operators
/// @param op
/// @param uval
/// @param umod
extern "C" void _ciface_quiccir_jw_prj_complexf64_DCCSC3D_complexf64_DCCSC3D(void* obj, view3_cd_t* pUval, view3_cd_t* pUmod)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_jw_prj_complexf64_DCCSC3D_complexf64_DCCSC3D\n";
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
    assert(pUmod->pos == pUval->pos);
    assert(pUmod->coo == pUval->coo);
    Tin viewMod(pUmod->data, pUmod->dataSize, pUmod->dims, pointers, indices);
    Tout viewVal(pUval->data, pUval->dataSize, pUval->dims, pointers, indices);
    // Check that op was set up
    auto cl = reinterpret_cast<op_t*>(obj);
    if (cl->getOp().data() == nullptr)
    {
        // dim 0 - Nr - radial points
        // dim 1 - N  - radial modes
        // dim 2 - L  - harmonic degree
        std::array<std::uint32_t, rank> dims {pUval->dims[0], pUmod->dims[0], pUmod->dims[2]};
        std::vector<std::uint32_t> layers;
        for (std::size_t i = 0; i < dims[2]; ++i) {
            if (pUmod->pos[i+1]-pUmod->pos[i] > 0)
            {
                layers.push_back(i);
            }
        }
        cl->allocOp(dims, layers);
        /// Set grid \todo set once per operator kind
        ::QuICC::Internal::Array igrid;
        ::QuICC::Internal::Array iweights;
        ::QuICC::Polynomial::Quadrature::WorlandRule quad;
        quad.computeQuadrature(igrid, iweights, pUval->dims[0]);
        // Populate op
        auto opView = cl->getOp();
        using namespace QuICC::Polynomial::Worland;
        using namespace QuICC::Transform::Worland;
        Builder<Top, QuICC::DenseSM::Worland::Operator<Wnl>, bwd_t> tBuilder;
        tBuilder.compute(opView, igrid, ::QuICC::Internal::Array());
    }
    assert(cl->getOp().data() != nullptr);
    // call
    cl->apply(viewVal, viewMod);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a jw prj operator
/// row major data, gpu operators
/// @param op
/// @param uval
/// @param umod
extern "C" void _ciface_quiccir_jw_prj_complexf64_DCCSC3DJIK_complexf64_DCCSC3DJIK(void* obj, view3_cd_t* pUval, view3_cd_t* pUmod)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_jw_prj_complexf64_DCCSC3DJIK_complexf64_DCCSC3DJIK\n";
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
    using namespace QuICC::Transform::Quadrature;
    using Tin = C_DCCSC3DJIK_t;
    using Tout = C_DCCSC3DJIK_t;
    using Top = QuICC::View::View<double, QuICC::View::CSL3D>;
    using backend_t = Cuda::ImplOp<Tout, Tin, Top>;
    using op_t = Op<Tout, Tin, Top, backend_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    ViewBase<std::uint32_t> pointers[rank];
    pointers[1] = ViewBase<std::uint32_t>(pUmod->pos, pUmod->posSize);
    ViewBase<std::uint32_t> indices[rank];
    indices[1] = ViewBase<std::uint32_t>(pUmod->coo, pUmod->cooSize);
    assert(pUmod->pos == pUval->pos);
    assert(pUmod->coo == pUval->coo);
    Tin viewMod(pUmod->data, pUmod->dataSize, pUmod->dims, pointers, indices);
    Tout viewVal(pUval->data, pUval->dataSize, pUval->dims, pointers, indices);
    // Check that op was set up
    auto cl = reinterpret_cast<op_t*>(obj);
    if (cl->getOp().data() == nullptr)
    {
        // dim 0 - Nr - radial points
        // dim 1 - N  - radial modes
        // dim 2 - L  - harmonic degree
        std::array<std::uint32_t, rank> dims {pUval->dims[0], pUmod->dims[0], pUmod->dims[2]};
        std::vector<std::uint32_t> layers;
        /// \todo getLayers
        for (std::size_t i = 0; i < dims[2]; ++i) {
            layers.push_back(i);
        }
        cl->allocOp(dims, layers);
        /// Set grid \todo set once per operator kind
        ::QuICC::Internal::Array igrid;
        ::QuICC::Internal::Array iweights;
        ::QuICC::Polynomial::Quadrature::WorlandRule quad;
        quad.computeQuadrature(igrid, iweights, pUval->dims[0]);
        // Populate op
        auto opView = cl->getOp();
        using namespace QuICC::Polynomial::Worland;
        using namespace QuICC::Transform::Worland;
        Builder<Top, QuICC::DenseSM::Worland::Operator<Wnl>, bwd_t> tBuilder;
        tBuilder.compute(opView, igrid, ::QuICC::Internal::Array());
    }
    assert(cl->getOp().data() != nullptr);
    // call
    cl->apply(viewVal, viewMod);
};
#endif
