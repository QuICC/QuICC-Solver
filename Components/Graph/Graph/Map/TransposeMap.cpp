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

void MapOps::setTranspose(mlir::quiccir::TransposeOp op)
{
    // Get index from MLIR source
    std::uint64_t index = op.getImplptr().value();
    if (index >= _thisArr.size()) {
        _thisArr.resize(index+1, nullptr);
    }
    if (_thisArr[index] == nullptr) {
        // check perm attribute
        auto perm = op.getPermutation();
        if (_isCpu)
        {
          using namespace QuICC::Transpose::Cpu;
          using namespace QuICC::Transpose;
          // check type attributes
          using namespace mlir;
          Type inTy = op.getInput().getType();
          auto inTensor = inTy.cast<RankedTensorType>();
          std::string inTyStr = inTensor.getEncoding().cast<StringAttr>().str();
          Type outTy = op.getOutput().getType();
          auto outTensor = outTy.cast<RankedTensorType>();
          std::string outTyStr = outTensor.getEncoding().cast<StringAttr>().str();
          if (outTyStr == "C_DCCSC3D_t" &&
            inTyStr == "C_DCCSC3D_t" &&
            perm[0] == 2 && perm[1] == 0)
          {
            using Tout = C_DCCSC3D_t;
            using Tin = C_DCCSC3D_t;
            using op_t = Op<Tout, Tin, p201_t>;
            _ops.push_back(std::make_unique<op_t>());
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
          }
          else if (outTyStr == "C_DCCSC3D_t" &&
            inTyStr == "C_DCCSC3D_t" &&
            perm[0] == 1 && perm[1] == 2)
          {
            using Tout = C_DCCSC3D_t;
            using Tin = C_DCCSC3D_t;
            using op_t = Op<Tout, Tin, p120_t>;
            _ops.push_back(std::make_unique<op_t>());
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
          }
          else if (outTyStr == "C_DCCSC3D_t" &&
            inTyStr == "C_S1CLCSC3D_t" &&
            perm[0] == 2 && perm[1] == 0)
          {
            using Tout = C_DCCSC3D_t;
            using Tin = C_S1CLCSC3D_t;
            using op_t = Op<Tout, Tin, p201_t>;
            _ops.push_back(std::make_unique<op_t>());
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
          }
          else if (outTyStr == "C_S1CLCSC3D_t" &&
            inTyStr == "C_DCCSC3D_t" &&
            perm[0] == 1 && perm[1] == 2)
          {
            using Tout = C_S1CLCSC3D_t;
            using Tin = C_DCCSC3D_t;
            using op_t = Op<Tout, Tin, p120_t>;
            _ops.push_back(std::make_unique<op_t>());
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
          }
        }
        #ifdef QUICC_HAS_CUDA_BACKEND
        else
        {
          using namespace QuICC::Transpose::Cuda;
          using namespace QuICC::Transpose;
          // check type attribute
          using namespace mlir;
          Type inTy = op.getInput().getType();
          auto inTensor = inTy.cast<RankedTensorType>();
          std::string inTyStr = inTensor.getEncoding().cast<StringAttr>().str();
          Type outTy = op.getOutput().getType();
          auto outTensor = outTy.cast<RankedTensorType>();
          std::string outTyStr = outTensor.getEncoding().cast<StringAttr>().str();
          if (outTyStr == "C_DCCSC3DJIK_t" &&
            inTyStr == "C_DCCSC3D_t" &&
            perm[0] == 2 && perm[1] == 0)
          {
            using Tout = C_DCCSC3DJIK_t;
            using Tin = C_DCCSC3D_t;
            using op_t = Op<Tout, Tin, p201_t>;
            _ops.push_back(std::make_unique<op_t>());
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
          }
          else if (outTyStr == "C_DCCSC3D_t" &&
            inTyStr == "C_DCCSC3DJIK_t" &&
            perm[0] == 1 && perm[2] == 0)
          {
            using Tout = C_DCCSC3D_t;
            using Tin = C_DCCSC3DJIK_t;
            using op_t = Op<Tout, Tin, p120_t>;
            _ops.push_back(std::make_unique<op_t>());
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
          }
          else if (outTyStr == "C_DCCSC3DJIK_t" &&
            inTyStr == "C_S1CLCSC3DJIK_t" &&
            perm[0] == 2 && perm[1] == 0)
          {
            using Tout = C_DCCSC3DJIK_t;
            using Tin = C_S1CLCSC3DJIK_t;
            using op_t = Op<Tout, Tin, p201_t>;
            _ops.push_back(std::make_unique<op_t>());
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
          }
          else if (outTyStr == "C_S1CLCSC3DJIK_t" &&
            inTyStr == "C_DCCSC3DJIK_t" &&
            perm[0] == 1 && perm[1] == 2)
          {
            using Tout = C_S1CLCSC3DJIK_t;
            using Tin = C_DCCSC3DJIK_t;
            using op_t = Op<Tout, Tin, p120_t>;
            _ops.push_back(std::make_unique<op_t>());
            auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
            assert(ptr != nullptr);
            _thisArr[index] = ptr;
          }
        }
        #endif
    }
    else {
        #ifndef NDEBUG
        std::cout << "transpose operator already allocated\n";
        #endif
    }
}

} // namespace Graph
} // namespace QuICC
