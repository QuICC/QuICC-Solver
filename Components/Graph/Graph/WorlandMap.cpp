// External includes
//
#include <iostream>

// Project includes
//
#include "Graph/OpsMap.hpp"
#include "Graph/BackendsMap.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandRule.hpp"
#include "ViewOps/Worland/Tags.hpp"
#include "ViewOps/Worland/Builder.hpp"
#include "DenseSM/Worland/Operator.hpp"


namespace QuICC
{
namespace Graph
{

void MapOps::setWorlandPrj(mlir::quiccir::JWPOp op)
{
    assert(false && "not implemented");
}

void MapOps::setWorlandInt(mlir::quiccir::JWIOp op)
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
            using Tout = C_DCCSC3D_t;
            using Top = QuICC::View::View<double, QuICC::View::CSL3DJIK>;
            using backend_t = Cpu::ImplOp<Tout, Tin, Top>;
            using op_t = Op<Tout, Tin, Top, backend_t>;
            _ops.push_back(std::make_unique<op_t>(_mem));
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            // Setup operator
            auto jwIntOp = dynamic_cast<op_t*>(ptr);
            constexpr size_t rank = 3;
            // dim 0 - N  - radial modes
            // dim 1 - Nr - radial points
            // dim 2 - L  - harmonic degree
            std::array<std::uint32_t, rank> dims {_modsDims[2], _physDims[2], _modsDims[1]};
            std::vector<std::uint32_t> layers;
            /// Dense operator \todo generalize for distributed op
            for (std::size_t i = 0; i < dims[2]; ++i) {
            layers.push_back(i);
            }
            jwIntOp->allocOp(dims, layers);
            /// Set grid \todo set once per operator kind
            Internal::Array igrid;
            Internal::Array iweights;
            ::QuICC::Polynomial::Quadrature::WorlandRule quad;
            quad.computeQuadrature(igrid, iweights, dims[1]);
            // Populate op
            auto opView = jwIntOp->getOp();
            using namespace QuICC::Polynomial::Worland;
            using namespace QuICC::Transform::Worland;
            Builder<Top, QuICC::DenseSM::Worland::Operator<Wnl>, fwd_t> tBuilder;
            tBuilder.compute(opView, igrid, iweights);
            // Add to thisArr
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
        }
        #ifdef QUICC_HAS_CUDA_BACKEND
        else
        {
            using namespace QuICC::Transform::Quadrature;
            using Tin = C_DCCSC3D_t;
            using Tout = C_DCCSC3D_t;
            using Top = QuICC::View::View<double, QuICC::View::CSL3D>;
            using backend_t = Cpu::ImplOp<Tout, Tin, Top>;
            using op_t = Op<Tout, Tin, Top, backend_t>;
            _ops.push_back(std::make_unique<op_t>(_mem));
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            // Setup operator
            auto jwIntOp = dynamic_cast<op_t*>(ptr);
            constexpr size_t rank = 3;
            // dim 0 - N  - radial modes
            // dim 1 - Nr - radial points
            // dim 2 - L  - harmonic degree
            std::array<std::uint32_t, rank> dims {_modsDims[2], _physDims[2], _modsDims[1]};
            std::vector<std::uint32_t> layers;
            /// Dense operator \todo generalize for distributed op
            for (std::size_t i = 0; i < dims[2]; ++i) {
            layers.push_back(i);
            }
            jwIntOp->allocOp(dims, layers);
            /// Set grid \todo set once per operator kind
            Internal::Array igrid;
            Internal::Array iweights;
            ::QuICC::Polynomial::Quadrature::WorlandRule quad;
            quad.computeQuadrature(igrid, iweights, dims[1]);
            // Populate op
            auto opView = jwIntOp->getOp();
            using namespace QuICC::Polynomial::Worland;
            using namespace QuICC::Transform::Worland;
            Builder<Top, QuICC::DenseSM::Worland::Operator<Wnl>, fwd_t> tBuilder;
            tBuilder.compute(opView, igrid, iweights);
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
