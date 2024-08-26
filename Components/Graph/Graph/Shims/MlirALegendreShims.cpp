#include <iostream>
#include <complex>
#include <cassert>

#include "Graph/BackendsMap.hpp"
#include "Graph/Types.hpp"
#include "Graph/Shims/MlirShims.hpp"

#include "QuICC/Polynomial/ALegendre/Plm.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
#include "ViewOps/ALegendre/Builder.hpp"
#include "Types/Internal/Math.hpp"


using namespace QuICC::Graph;

/// @brief C Interface to MLIR for a al int operator
/// column major data, cpu operators
/// @param op
/// @param uval
/// @param umod
extern "C" void _ciface_quiccir_al_int_complexf64_S1CLCSC3D_complexf64_DCCSC3D(void* obj, view3_cd_t* pUmod, view3_cd_t* pUval)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_al_int_complexf64_S1CLCSC3D_complexf64_DCCSC3D\n";
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
    assert(pUmod->pos == pUval->pos);
    assert(pUmod->coo == pUval->coo);
    Tin viewVal(pUval->data, pUval->dataSize, pUval->dims, pointers, indices);
    Tout viewMod(pUmod->data, pUmod->dataSize, pUmod->dims, pointers, indices);
    // Check that op was set up
    auto cl = reinterpret_cast<op_t*>(obj);
    if (cl->getOp().data() == nullptr)
    {
        /// dim 0 - L  - harmonic degree
        /// dim 1 - Ntheta - colatitudinal points
        /// dim 2 - M  - harmonic order
        std::array<std::uint32_t, rank> dims {pUmod->dims[0], pUval->dims[0], pUmod->dims[2]};
        std::vector<std::uint32_t> layers;
        /// Dense operator \todo generalize for distributed op
        for (std::size_t i = 0; i < dims[2]; ++i) {
            layers.push_back(i);
        }
        cl->allocOp(dims, layers);
        /// Set grid \todo set once per operator kind
        ::QuICC::Internal::Array igrid;
        ::QuICC::Internal::Array iweights;
        ::QuICC::Polynomial::Quadrature::LegendreRule quad;
        quad.computeQuadrature(igrid, iweights, pUval->dims[0]);
        // scale for spherical harmonics
        iweights.array() *= 2.0*::QuICC::Internal::Math::PI;
        // Populate op
        auto opView = cl->getOp();
        using namespace QuICC::Transform::ALegendre;
        builder<Top, ::QuICC::Polynomial::ALegendre::Plm, ::QuICC::Internal::Array::Scalar, 0>(opView, igrid, iweights);
    }
    assert(cl->getOp().data() != nullptr);
    // call
    cl->apply(viewMod, viewVal);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a al int operator
/// row major, gpu operators
/// @param op
/// @param uval
/// @param umod
extern "C" void _ciface_quiccir_al_int_C_S1CLCSC3DJIK_t_C_DCCSC3DJIK_t(void* obj, view3_cd_t* pUmod, view3_cd_t* pUval)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_al_int_C_S1CLCSC3DJIK_t_C_DCCSC3DJIK_t\n";
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
    using Tout = C_S1CLCSC3DJIK_t;
    using Top = QuICC::View::View<double, QuICC::View::S1CLCSC3D>;
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
    Tin viewVal(pUval->data, pUval->dataSize, pUval->dims, pointers, indices);
    Tout viewMod(pUmod->data, pUmod->dataSize, pUmod->dims, pointers, indices);
    // Check that op was set up
    auto cl = reinterpret_cast<op_t*>(obj);
    if (cl->getOp().data() == nullptr)
    {
        /// dim 0 - L  - harmonic degree
        /// dim 1 - Ntheta - colatitudinal points
        /// dim 2 - M  - harmonic order
        std::array<std::uint32_t, rank> dims {pUmod->dims[0], pUval->dims[0], pUmod->dims[2]};
        std::vector<std::uint32_t> layers;
        /// Dense operator \todo generalize for distributed op
        for (std::size_t i = 0; i < dims[2]; ++i) {
            layers.push_back(i);
        }
        cl->allocOp(dims, layers);
        /// Set grid \todo set once per operator kind
        ::QuICC::Internal::Array igrid;
        ::QuICC::Internal::Array iweights;
        ::QuICC::Polynomial::Quadrature::LegendreRule quad;
        quad.computeQuadrature(igrid, iweights, pUval->dims[0]);
        // scale for spherical harmonics
        iweights.array() *= 2.0*::QuICC::Internal::Math::PI;
        // Populate op
        auto opView = cl->getOp();
        using namespace QuICC::Transform::ALegendre;
        builder<Top, ::QuICC::Polynomial::ALegendre::Plm, ::QuICC::Internal::Array::Scalar, 0>(opView, igrid, iweights);
    }
    assert(cl->getOp().data() != nullptr);
    // call
    cl->apply(viewMod, viewVal);
};
#endif

/// @brief C Interface to MLIR for a al prj operator
/// column major, cpu operators
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_al_prj_complexf64_DCCSC3D_complexf64_S1CLCSC3D(void* obj,  view3_cd_t* pUval, view3_cd_t* pUmod)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_al_prj_complexf64_DCCSC3D_complexf64_S1CLCSC3D\n";
    #endif
    assert(obj != nullptr);
    assert(pUval != nullptr);
    assert(pUmod != nullptr);
    // op
    using namespace QuICC::Transform::Quadrature;
    using Tin = C_S1CLCSC3D_t;
    using Tout = C_DCCSC3D_t;
    using Top = QuICC::View::View<double, QuICC::View::CS1RL3DJIK>;
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
        /// dim 0 - Ntheta - colatitudinal points
        /// dim 1 - L  - harmonic degree
        /// dim 2 - M  - harmonic order
        std::array<std::uint32_t, rank> dims {pUval->dims[0], pUmod->dims[0], pUmod->dims[2]};
        std::vector<std::uint32_t> layers;
        /// Dense operator \todo generalize for distributed op
        for (std::size_t i = 0; i < dims[2]; ++i) {
            layers.push_back(i);
        }
        cl->allocOp(dims, layers);
        /// Set grid \todo set once per operator kind
        ::QuICC::Internal::Array igrid;
        ::QuICC::Internal::Array iweights;
        ::QuICC::Polynomial::Quadrature::LegendreRule quad;
        quad.computeQuadrature(igrid, iweights, pUval->dims[0]);
        // Populate op
        auto opView = cl->getOp();
        using namespace QuICC::Transform::ALegendre;
        builder<Top, ::QuICC::Polynomial::ALegendre::Plm, ::QuICC::Internal::Array::Scalar, 0>(opView, igrid, ::QuICC::Internal::Array());
    }
    assert(cl->getOp().data() != nullptr);
    // call
    cl->apply(viewVal, viewMod);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a al int operator
/// row major, gpu operators
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_al_prj_complexf64_DCCSC3DJIK_complexf64_S1CLCSC3DJIK(void* obj,  view3_cd_t* pUval, view3_cd_t* pUmod)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_al_prj_complexf64_DCCSC3DJIK_complexf64_S1CLCSC3DJIK\n";
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
    using Tin = C_S1CLCSC3DJIK_t;
    using Tout = C_DCCSC3DJIK_t;
    using Top = QuICC::View::View<double, QuICC::View::CS1RL3D>;
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
        /// dim 0 - Ntheta - colatitudinal points
        /// dim 1 - L  - harmonic degree
        /// dim 2 - M  - harmonic order
        std::array<std::uint32_t, rank> dims {pUval->dims[0], pUmod->dims[0], pUmod->dims[2]};
        std::vector<std::uint32_t> layers;
        /// Dense operator \todo generalize for distributed op
        for (std::size_t i = 0; i < dims[2]; ++i) {
            layers.push_back(i);
        }
        cl->allocOp(dims, layers);
        /// Set grid \todo set once per operator kind
        ::QuICC::Internal::Array igrid;
        ::QuICC::Internal::Array iweights;
        ::QuICC::Polynomial::Quadrature::LegendreRule quad;
        quad.computeQuadrature(igrid, iweights, pUval->dims[0]);
        // Populate op
        auto opView = cl->getOp();
        using namespace QuICC::Transform::ALegendre;
        builder<Top, ::QuICC::Polynomial::ALegendre::Plm, ::QuICC::Internal::Array::Scalar, 0>(opView, igrid, ::QuICC::Internal::Array());
    }
    assert(cl->getOp().data() != nullptr);
    // call
    cl->apply(viewVal, viewMod);
};
#endif
