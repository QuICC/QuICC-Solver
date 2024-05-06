// External includes
//
#include <iostream>

// Project includes
//
#include "Graph/OpsMap.hpp"
#include "Graph/BackendsMap.hpp"


namespace QuICC
{
namespace Graph
{

void MapOps::setFourierPrj(mlir::quiccir::FrPOp op)
{
    // Get index from MLIR source
    std::uint64_t index = op.getImplptr().value();
    if (index >= _thisArr.size()) {
        _thisArr.resize(index+1, nullptr);
    }
    if (_thisArr[index] == nullptr) {
        if (_isCpu)
        {
            using namespace QuICC::Transform::Fourier;
            using backend_t = QuICC::Graph::viewCpu_t;
            using Tin = C_DCCSC3D_t;
            using Tout = R_DCCSC3D_t;
            using backendFft_t = Fft_t<backend_t, Tout, Tin>;
            using backendDiff_t = MixedDiff_t<backend_t, Tin, 0, bwd_t,
                QuICC::Transform::Fourier::none_m>;
            using op_t = Mixed::Projector::DOp<Tout, Tin, backendFft_t,
            backendDiff_t>;
            _ops.push_back(std::make_unique<op_t>(_mem));
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
        }
        #ifdef QUICC_HAS_CUDA_BACKEND
        else
        {
            using namespace QuICC::Transform::Fourier;
            using backend_t = QuICC::Graph::viewGpu_t;
            using Tin = C_DCCSC3D_t;
            using Tout = R_DCCSC3D_t;
            using backendFft_t = Fft_t<backend_t, Tout, Tin>;
            using backendDiff_t = MixedDiff_t<backend_t, Tin, 0, bwd_t,
                QuICC::Transform::Fourier::none_m>;
            using op_t = Mixed::Projector::DOp<Tout, Tin, backendFft_t,
            backendDiff_t>;
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
    // Get index from MLIR source
    std::uint64_t index = op.getImplptr().value();
    if (index >= _thisArr.size()) {
        _thisArr.resize(index+1, nullptr);
    }
    if (_thisArr[index] == nullptr) {
        if (_isCpu)
        {
            using namespace QuICC::Transform::Fourier;
            using backend_t = QuICC::Graph::viewCpu_t;
            using Tin = R_DCCSC3D_t;
            using Tout = C_DCCSC3D_t;
            using backendFft_t = Fft_t<backend_t, Tout, Tin>;
            using backendDiff_t = MixedDiff_t<backend_t, Tout, 0, fwd_t,
                QuICC::Transform::Fourier::none_m>;
            using op_t = Mixed::Integrator::DOp<Tout, Tin, backendFft_t,
            backendDiff_t>;
            _ops.push_back(std::make_unique<op_t>());
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
        }
        #ifdef QUICC_HAS_CUDA_BACKEND
        else
        {
            using namespace QuICC::Transform::Fourier;
            using backend_t = QuICC::Graph::viewGpu_t;
            using Tin = R_DCCSC3D_t;
            using Tout = C_DCCSC3D_t;
            using backendFft_t = Fft_t<backend_t, Tout, Tin>;
            using backendDiff_t = MixedDiff_t<backend_t, Tout, 0, fwd_t,
                QuICC::Transform::Fourier::none_m>;
            using op_t = Mixed::Integrator::DOp<Tout, Tin, backendFft_t,
            backendDiff_t>;
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
