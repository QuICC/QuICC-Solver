// External includes
//
#include <iostream>

// Project includes
//
#include "Graph/OpsMap.hpp"
#include "ViewOps/Fourier/Mixed/OpsType.hpp"

namespace QuICC
{
namespace Graph
{

void MapOps::setFourierPrj(mlir::quiccir::FrPOp op)
{
    using namespace QuICC::Transform::Fourier;
    // Get index from MLIR source
    std::uint64_t index = op.getImplptr().value();
    if (index >= _thisArr.size()) {
        _thisArr.resize(index+1, nullptr);
    }
    if (_thisArr[index] == nullptr) {
        using Tin = C_DCCSC3D_t;
        using Tout = R_DCCSC3D_t;
        /// \todo check kind
        if (_isCpu)
        {
            using backend_t = viewCpu_t;
            using op_t = Mixed::OpsType<Tout, Tin, P_t, bwd_t, backend_t>;
            _ops.push_back(std::make_unique<op_t>(_mem));
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
        }
        #ifdef QUICC_HAS_CUDA_BACKEND
        else
        {
            using backend_t = viewGpu_t;
            using op_t = Mixed::OpsType<Tout, Tin, P_t, bwd_t, backend_t>;
            _ops.push_back(std::make_unique<op_t>(_mem));
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
        }
        #endif
    }
    else {
        #ifndef NDEBUG
        std::cout << "operator already allocated\n";
        #endif
    }
}

void MapOps::setFourierInt(mlir::quiccir::FrIOp op)
{
    using namespace QuICC::Transform::Fourier;
    // Get index from MLIR source
    std::uint64_t index = op.getImplptr().value();
    if (index >= _thisArr.size()) {
        _thisArr.resize(index+1, nullptr);
    }
    if (_thisArr[index] == nullptr) {
        using Tin = R_DCCSC3D_t;
        using Tout = C_DCCSC3D_t;
        if (_isCpu)
        {
            using backend_t = viewCpu_t;
            //// \todo map each operator or split FFT and diff
            using op_t = Mixed::OpsType<Tout, Tin, P_t, fwd_t, backend_t>;
            _ops.push_back(std::make_unique<op_t>());
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
        }
        #ifdef QUICC_HAS_CUDA_BACKEND
        else
        {
            using backend_t = viewGpu_t;
            using op_t = Mixed::OpsType<Tout, Tin, P_t, fwd_t, backend_t>;
            _ops.push_back(std::make_unique<op_t>());
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
        }
        #endif
    }
    else {
        #ifndef NDEBUG
        std::cout << "operator already allocated\n";
        #endif
    }
}

} // namespace Graph
} // namespace QuICC
