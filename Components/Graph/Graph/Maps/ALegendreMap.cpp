// External includes
//
#include <iostream>

// Project includes
//
#include "Graph/OpsMap.hpp"
#include "ViewOps/Quadrature/Op.hpp"
#include "ViewOps/Quadrature/Impl.hpp"

namespace QuICC
{
namespace Graph
{

void MapOps::setALegendrePrj(mlir::quiccir::AlPOp op)
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
            using Tin = C_S1CLCSC3D_t;
            using Tout = C_DCCSC3D_t;
            using Top = QuICC::View::View<double, QuICC::View::CS1RL3DJIK>;
            if (op.getKind() == "DivS1Dp" ||
                op.getKind() == "LlDivS1Dp")
            {
                // Special treatment operators
                // aka fused Phi diff
                using backend_t = Cpu::ImplOp<Tout, Tin, Top, Transform::Quadrature::diffPhiPrj_m>;
                using op_t = Op<Tout, Tin, Top, backend_t>;
                _ops.push_back(std::make_unique<op_t>(_mem));
                auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
                // Add to thisArr
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            else
            {
                // Regular AL ops
                using backend_t = Cpu::ImplOp<Tout, Tin, Top, Transform::Quadrature::none_m>;
                using op_t = Op<Tout, Tin, Top, backend_t>;
                _ops.push_back(std::make_unique<op_t>(_mem));
                auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
                // Add to thisArr
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
        }
        #ifdef QUICC_HAS_CUDA_BACKEND
        else
        {
            using namespace QuICC::Transform::Quadrature;
            using Tin = C_S1CLCSC3DJIK_t;
            using Tout = C_DCCSC3DJIK_t;
            using Top = QuICC::View::View<double, QuICC::View::CS1RL3D>;
            using QuICC::Transform::Quadrature::Cuda::ImplOp;
            if (op.getKind() == "DivS1Dp" ||
                op.getKind() == "LlDivS1Dp")
            {
                // Special treatment operators
                // aka fused Phi diff
                using backend_t = ImplOp<Tout, Tin, Top, Transform::Quadrature::diffPhiPrj_m>;
                using op_t = Op<Tout, Tin, Top, backend_t>;
                _ops.push_back(std::make_unique<op_t>(_mem));
                auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
                // Add to thisArr
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            else
            {
                // Regular AL Ops
                using backend_t = ImplOp<Tout, Tin, Top, Transform::Quadrature::none_m>;
                using op_t = Op<Tout, Tin, Top, backend_t>;
                _ops.push_back(std::make_unique<op_t>(_mem));
                auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
                // Add to thisArr
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
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
            if (op.getKind() == "DivS1Dp" ||
                op.getKind() == "LlDivS1Dp"||
                op.getKind() == "DivLlDivS1Dp")
            {
                // Special treatment operators
                // aka fused Phi diff
                using backend_t = Cpu::ImplOp<Tout, Tin, Top, Transform::Quadrature::diffPhiInt_m>;
                using op_t = Op<Tout, Tin, Top, backend_t>;
                _ops.push_back(std::make_unique<op_t>(_mem));
                auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
                // Add to thisArr
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            else
            {
                // Regular AL Ops
                using backend_t = Cpu::ImplOp<Tout, Tin, Top>;
                using op_t = Op<Tout, Tin, Top, backend_t>;
                _ops.push_back(std::make_unique<op_t>(_mem));
                auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
                // Add to thisArr
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }


        }
        #ifdef QUICC_HAS_CUDA_BACKEND
        else
        {
            using namespace QuICC::Transform::Quadrature;
            using Tin = C_DCCSC3DJIK_t;
            using Tout = C_S1CLCSC3DJIK_t;
            using Top = QuICC::View::View<double, QuICC::View::S1CLCSC3D>;
            using QuICC::Transform::Quadrature::Cuda::ImplOp;
            if (op.getKind() == "DivS1Dp" ||
                op.getKind() == "LlDivS1Dp"||
                op.getKind() == "DivLlDivS1Dp")
            {
                // Special treatment operators
                // aka fused Phi diff
                using backend_t = ImplOp<Tout, Tin, Top, Transform::Quadrature::diffPhiInt_m>;
                using op_t = Op<Tout, Tin, Top, backend_t>;
                _ops.push_back(std::make_unique<op_t>(_mem));
                auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
                // Add to thisArr
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            else
            {
                // Regular AL Ops
                using backend_t = ImplOp<Tout, Tin, Top>;
                using op_t = Op<Tout, Tin, Top, backend_t>;
                _ops.push_back(std::make_unique<op_t>(_mem));
                auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
                // Add to thisArr
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
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
