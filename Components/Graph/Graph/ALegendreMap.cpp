// External includes
//

// Project includes
//
#include "Graph/OpsMap.hpp"
#include "Graph/BackendsMap.hpp"
#include "QuICC/Polynomial/ALegendre/Plm.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
#include "ViewOps/ALegendre/Builder.hpp"


namespace QuICC
{
namespace Graph
{

void MapOps::setALegendrePrj(mlir::quiccir::AlPOp op)
{
    assert(false && "not implemented");
}

void MapOps::setALegendreInt(mlir::quiccir::AlIOp op)
{
    // Get index from MLIR source
    std::uint64_t index = op.getImplptr().value();
    if (index >= _thisArr.size())
    {
        _thisArr.resize(index+1, nullptr);
    }
    if (_thisArr[index] == nullptr)
    {
        if (_isCpu)
        {
            using namespace QuICC::Transform::Quadrature;
            using Tin = C_DCCSC3D_t;
            using Tout = C_S1CLCSC3D_t;
            using Top = QuICC::View::View<double, QuICC::View::S1CLCSC3DJIK>;
            using backend_t = Cpu::ImplOp<Tout, Tin, Top>;
            using op_t = Op<Tout, Tin, Top, backend_t>;
            _ops.push_back(std::make_unique<op_t>(_mem));
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            // Setup operator
            auto alIntOp = dynamic_cast<op_t*>(ptr);
            constexpr size_t rank = 3;
            /// dim 0 - L  - harmonic degree
            /// dim 1 - Ntheta - colatitudinal points
            /// dim 2 - M  - harmonic order
            std::array<std::uint32_t, rank> dims {_modsDims[1], _physDims[1], _modsDims[0]};
            std::vector<std::uint32_t> layers;
            /// Dense operator \todo generalize for distributed op
            for (std::size_t i = 0; i < dims[2]; ++i) {
                layers.push_back(i);
            }
            alIntOp->allocOp(dims, layers);
            /// Set grid \todo set once per operator kind
            Internal::Array igrid;
            Internal::Array iweights;
            ::QuICC::Polynomial::Quadrature::LegendreRule quad;
            quad.computeQuadrature(igrid, iweights, dims[1]);
            // Populate op
            auto opView = alIntOp->getOp();
            using namespace QuICC::Transform::ALegendre;
            builder<Top, Polynomial::ALegendre::Plm, Internal::Array::Scalar, 0>(opView, igrid, iweights);
            // Add to thisArr
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
        }
        #ifdef QUICC_HAS_CUDA_BACKEND
        else
        {
            using namespace QuICC::Transform::Quadrature;
            using Tin = C_DCCSC3DJIK_t;
            using Tout = C_S1CLCSC3DJIK_t;
            using Top = QuICC::View::View<double, QuICC::View::S1CLCSC3D>;
            using backend_t = Cuda::ImplOp<Tout, Tin, Top>;
            using op_t = Op<Tout, Tin, Top, backend_t>;
            _ops.push_back(std::make_unique<op_t>(_mem));
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            // Setup operator
            auto alIntOp = dynamic_cast<op_t*>(ptr);
            constexpr size_t rank = 3;
            /// dim 0 - L  - harmonic degree
            /// dim 1 - Ntheta - colatitudinal points
            /// dim 2 - M  - harmonic order
            std::array<std::uint32_t, rank> dims {_modsDims[1], _physDims[1], _modsDims[0]};
            std::vector<std::uint32_t> layers;
            /// Dense operator \todo generalize for distributed op
            for (std::size_t i = 0; i < dims[2]; ++i) {
                layers.push_back(i);
            }
            alIntOp->allocOp(dims, layers);
            /// Set grid \todo set once per operator kind
            Internal::Array igrid;
            Internal::Array iweights;
            ::QuICC::Polynomial::Quadrature::LegendreRule quad;
            quad.computeQuadrature(igrid, iweights, dims[1]);
            // Populate op
            auto opView = alIntOp->getOp();
            using namespace QuICC::Transform::ALegendre;
            builder<Top, Polynomial::ALegendre::Plm, Internal::Array::Scalar, 0>(opView, igrid, iweights);
            // Add to thisArr
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
        }
        #endif
    }
    else
    {
        #ifndef NDEBUG
        std::cout << "operator already allocated\n";
        #endif
    }
}

} // namespace Graph
} // namespace QuICC
