// External includes
//
#include <iostream>

// Project includes
//
#include "Graph/OpsMap.hpp"
#include "ViewOps/Slicewise/Op.hpp"
#include "ViewOps/Slicewise/Functors.hpp"

#include "QuICC/Polynomial/Quadrature/WorlandRule.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"

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
        using T = R_DCCSC3D_t;
        using namespace QuICC::Slicewise;
        double scaling = 1.0;
        if (op.getKind() == "coriolis_cos")
        {
            scaling = _physParams.coriolis;
            /// \todo check type
            if (_isCpu)
            {
                using namespace QuICC::Slicewise::Cpu;
                using op_t = Op<1, QuICC::Polynomial::Quadrature::LegendreRule, MulCosFunctor<double>, T, T>;
                _ops.push_back(std::make_unique<op_t>(MulCosFunctor<double>(scaling), _mem));
                auto* ptr = std::get<std::shared_ptr<NaryOp<T, T>>>(_ops.back()).get();
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            #ifdef QUICC_HAS_CUDA_BACKEND
            else
            {
                using namespace QuICC::Slicewise::Cuda;
                using op_t = Op<1, QuICC::Polynomial::Quadrature::LegendreRule, MulCosFunctor<double>, T, T>;
                _ops.push_back(std::make_unique<op_t>(MulCosFunctor<double>(scaling), _mem));
                auto* ptr = std::get<std::shared_ptr<NaryOp<T, T>>>(_ops.back()).get();
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            #endif
        }
        else if (op.getKind() == "coriolis_sin")
        {
            scaling = _physParams.coriolis;
            /// \todo check type
            if (_isCpu)
            {
                using namespace QuICC::Slicewise::Cpu;
                using op_t = Op<1, QuICC::Polynomial::Quadrature::LegendreRule, MulSinFunctor<double>, T, T>;
                _ops.push_back(std::make_unique<op_t>(MulSinFunctor<double>(scaling), _mem));
                auto* ptr = std::get<std::shared_ptr<NaryOp<T, T>>>(_ops.back()).get();
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            #ifdef QUICC_HAS_CUDA_BACKEND
            else
            {
                using namespace QuICC::Slicewise::Cuda;
                using op_t = Op<1, QuICC::Polynomial::Quadrature::LegendreRule, MulSinFunctor<double>, T, T>;
                _ops.push_back(std::make_unique<op_t>(MulSinFunctor<double>(scaling), _mem));
                auto* ptr = std::get<std::shared_ptr<NaryOp<T, T>>>(_ops.back()).get();
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            #endif
        }
        else
        {
            /// \todo check all kinds
            if (op.getKind() == "buoyancy")
            {
                scaling = _physParams.buoyancy;
            }
            else if (op.getKind() == "inertia")
            {
                scaling = _physParams.inertia;
            }
            else if (op.getKind() == "transport")
            {
                scaling = _physParams.transport;
            }
            /// \todo check type
            if (_isCpu)
            {
                using namespace QuICC::Slicewise::Cpu;
                using op_t = Op<2, QuICC::Polynomial::Quadrature::WorlandRule, MulRFunctor<double>, T, T>;
                _ops.push_back(std::make_unique<op_t>(MulRFunctor<double>(scaling), _mem));
                auto* ptr = std::get<std::shared_ptr<NaryOp<T, T>>>(_ops.back()).get();
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            #ifdef QUICC_HAS_CUDA_BACKEND
            else
            {
                using namespace QuICC::Slicewise::Cuda;
                using op_t = Op<2, QuICC::Polynomial::Quadrature::WorlandRule, MulRFunctor<double>, T, T>;
                _ops.push_back(std::make_unique<op_t>(MulRFunctor<double>(scaling), _mem));
                auto* ptr = std::get<std::shared_ptr<NaryOp<T, T>>>(_ops.back()).get();
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            #endif
        }
    }
    else {
        #ifndef NDEBUG
        std::cout << "operator already allocated\n";
        #endif
    }
}

} // namespace Graph
} // namespace QuICC
