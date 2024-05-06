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

void MapOps::setAdd(mlir::quiccir::AddOp op)
{
    // Get index from MLIR source
    std::uint64_t index = op.getImplptr().value();
    if (index >= _thisArr.size()) {
        _thisArr.resize(index+1, nullptr);
    }
    if (_thisArr[index] == nullptr) {
        if (_isCpu)
        {
            using namespace QuICC::Pointwise::Cpu;
            using namespace QuICC::Pointwise;
            using T = C_DCCSC3D_t;
            using op_t = Op<AddFunctor<std::complex<double>>, T, T, T>;
            _ops.push_back(std::make_unique<op_t>(AddFunctor<std::complex<double>>()));
            auto* ptr = std::get<std::shared_ptr<NaryOp<T, T, T>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
        }
        #ifdef QUICC_HAS_CUDA_BACKEND
        else
        {
            using namespace QuICC::Pointwise::Cuda;
            using namespace QuICC::Pointwise;
            using T = C_DCCSC3DJIK_t;
            using op_t = Op<AddFunctor<std::complex<double>>, T, T, T>;
            _ops.push_back(std::make_unique<op_t>(AddFunctor<std::complex<double>>()));
            auto* ptr = std::get<std::shared_ptr<NaryOp<T, T, T>>>(_ops.back()).get();
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

void MapOps::setSub(mlir::quiccir::SubOp op)
{
    // Get index from MLIR source
    std::uint64_t index = op.getImplptr().value();
    if (index >= _thisArr.size()) {
        _thisArr.resize(index+1, nullptr);
    }
    if (_thisArr[index] == nullptr) {
        if (_isCpu)
        {
            using namespace QuICC::Pointwise::Cpu;
            using namespace QuICC::Pointwise;
            using T = C_DCCSC3D_t;
            using op_t = Op<SubFunctor<std::complex<double>>, T, T, T>;
            _ops.push_back(std::make_unique<op_t>(SubFunctor<std::complex<double>>()));
            auto* ptr = std::get<std::shared_ptr<NaryOp<T, T, T>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
        }
        #ifdef QUICC_HAS_CUDA_BACKEND
        else
        {
            using namespace QuICC::Pointwise::Cuda;
            using namespace QuICC::Pointwise;
            using T = C_DCCSC3DJIK_t;
            using op_t = Op<SubFunctor<std::complex<double>>, T, T, T>;
            _ops.push_back(std::make_unique<op_t>(SubFunctor<std::complex<double>>()));
            auto* ptr = std::get<std::shared_ptr<NaryOp<T, T, T>>>(_ops.back()).get();
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
