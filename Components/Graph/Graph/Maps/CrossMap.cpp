// External includes
//
#include <iostream>

// Project includes
//
#include "Graph/OpsMap.hpp"
#include "ViewOps/Pointwise/Pointwise.hpp"
#include "ViewOps/Pointwise/Functors.hpp"


namespace QuICC
{
namespace Graph
{

void MapOps::setCross(mlir::quiccir::CrossOp op)
{
    // Get index from MLIR source
    std::uint64_t index = op.getImplptr().value();
    if (index >= _thisArr.size()) {
        _thisArr.resize(index+1, nullptr);
    }
    double scaling = 1.0;
    /// \todo check all kinds
    if (op.getKind() == "inertia")
    {
        scaling = _physParams.inertia;
    }
    if (_thisArr[index] == nullptr) {
        using namespace QuICC::Pointwise;
        using T = R_DCCSC3D_t;
        /// \todo check kind
        if (_isCpu)
        {
            using namespace QuICC::Pointwise::Cpu;
            using op_t = Op<CrossCompFunctor<double>, T, T, T, T, T>;
            _ops.push_back(std::make_unique<op_t>(CrossCompFunctor<double>(scaling)));
            auto* ptr = std::get<std::shared_ptr<NaryOp<T, T, T, T, T>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
        }
        #ifdef QUICC_HAS_CUDA_BACKEND
        else
        {
            using namespace QuICC::Pointwise::Cuda;
            using op_t = Op<CrossCompFunctor<double>, T, T, T, T, T>;
            _ops.push_back(std::make_unique<op_t>(CrossCompFunctor<double>(scaling)));
            auto* ptr = std::get<std::shared_ptr<NaryOp<T, T, T, T, T>>>(_ops.back()).get();
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
