// External includes
//
#include <iostream>

// Project includes
//
#include "Graph/OpsMap.hpp"
#include "ViewOps/Slicewise/Op.hpp"
#include "ViewOps/Slicewise/Functors.hpp"


namespace QuICC
{
namespace Graph
{

void MapOps::setMulConst(mlir::quiccir::MulConstOp op)
{
    // Get index from MLIR source
    std::uint64_t index = op.getImplptr().value();
    if (index >= _thisArr.size()) {
        _thisArr.resize(index+1, nullptr);
    }
    if (_thisArr[index] == nullptr) {
        /// \todo check all kinds
        /// \todo check type
        if (_isCpu)
        {
            double scaling = 1.0;
            if (op.getKind() == "buoyancy")
            {
                /// \todo set by Jitter
                scaling = 3122;
            }
            using namespace QuICC::Slicewise::Cpu;
            using namespace QuICC::Slicewise;
            using T = R_DCCSC3D_t;
            using op_t = Op<MulRFunctor<double>, T, T>;
            _ops.push_back(std::make_unique<op_t>(_mem, MulRFunctor<double>(scaling)));
            auto* ptr = std::get<std::shared_ptr<NaryOp<T, T>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
        }
        #ifdef QUICC_HAS_CUDA_BACKEND
        else
        {
            using namespace QuICC::Slicewise::Cuda;
            using namespace QuICC::Slicewise;
            using T = R_DCCSC3D_t;
            using op_t = Op<MulRFunctor<double>, T, T>;
            _ops.push_back(std::make_unique<op_t>(MulConstCompFunctor<double>()));
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
