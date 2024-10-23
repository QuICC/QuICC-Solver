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

void MapOps::setAdd(mlir::quiccir::AddOp op)
{
    // Get index from MLIR source
    /// \todo check encoding
    std::uint64_t index = op.getImplptr().value();
    if (index >= _thisArr.size()) {
        _thisArr.resize(index+1, nullptr);
    }
    if (_thisArr[index] == nullptr) {
        auto eleLhsTy = op.getLhs().getType().getElementType();
        auto eleRhsTy = op.getRhs().getType().getElementType();
        if (_isCpu)
        {
            using namespace QuICC::Pointwise::Cpu;
            using namespace QuICC::Pointwise;
            if (llvm::isa<mlir::Float64Type>(eleLhsTy) &&
                llvm::isa<mlir::Float64Type>(eleRhsTy))
            {
                using ptTy = double;
                using T = View::View<ptTy, View::DCCSC3D>;
                using op_t = Op<AddFunctor<ptTy>, T, T, T>;
                _ops.push_back(std::make_unique<op_t>(AddFunctor<ptTy>()));
                auto* ptr = std::get<std::shared_ptr<NaryOp<T, T, T>>>(_ops.back()).get();
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            else if (llvm::isa<mlir::ComplexType>(eleLhsTy) &&
                llvm::isa<mlir::Float64Type>(eleLhsTy.cast<mlir::ComplexType>().getElementType()) &&
                llvm::isa<mlir::ComplexType>(eleRhsTy) &&
                llvm::isa<mlir::Float64Type>(eleRhsTy.cast<mlir::ComplexType>().getElementType()))
            {
                using ptTy = std::complex<double>;
                using T = View::View<ptTy, View::DCCSC3D>;
                using op_t = Op<AddFunctor<ptTy>, T, T, T>;
                _ops.push_back(std::make_unique<op_t>(AddFunctor<ptTy>()));
                auto* ptr = std::get<std::shared_ptr<NaryOp<T, T, T>>>(_ops.back()).get();
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            else
            {
                throw std::logic_error("sub op not implemented.");
            }
        }
        #ifdef QUICC_HAS_CUDA_BACKEND
        else
        {
            using namespace QuICC::Pointwise::Cuda;
            using namespace QuICC::Pointwise;
            if (llvm::isa<mlir::Float64Type>(eleLhsTy) &&
                llvm::isa<mlir::Float64Type>(eleRhsTy))
            {
                using ptTy = double;
                using T = R_DCCSC3D_t;
                using op_t = Op<AddFunctor<ptTy>, T, T, T>;
                _ops.push_back(std::make_unique<op_t>(AddFunctor<ptTy>()));
                auto* ptr = std::get<std::shared_ptr<NaryOp<T, T, T>>>(_ops.back()).get();
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            else if (llvm::isa<mlir::ComplexType>(eleLhsTy) &&
                llvm::isa<mlir::Float64Type>(eleLhsTy.cast<mlir::ComplexType>().getElementType()) &&
                llvm::isa<mlir::ComplexType>(eleRhsTy) &&
                llvm::isa<mlir::Float64Type>(eleRhsTy.cast<mlir::ComplexType>().getElementType()))
            {
                using ptTy = std::complex<double>;
                using T = View::View<ptTy, View::DCCSC3D>;
                using op_t = Op<AddFunctor<ptTy>, T, T, T>;
                _ops.push_back(std::make_unique<op_t>(AddFunctor<ptTy>()));
                auto* ptr = std::get<std::shared_ptr<NaryOp<T, T, T>>>(_ops.back()).get();
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            else
            {
                throw std::logic_error("sub op not implemented.");
            }

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
        auto eleLhsTy = op.getLhs().getType().getElementType();
        auto eleRhsTy = op.getRhs().getType().getElementType();
        if (_isCpu)
        {
            using namespace QuICC::Pointwise::Cpu;
            using namespace QuICC::Pointwise;
            if (llvm::isa<mlir::Float64Type>(eleLhsTy) &&
                llvm::isa<mlir::Float64Type>(eleRhsTy))
            {
                using ptTy = double;
                using T = View::View<ptTy, View::DCCSC3D>;
                using op_t = Op<SubFunctor<ptTy>, T, T, T>;
                _ops.push_back(std::make_unique<op_t>(SubFunctor<ptTy>()));
                auto* ptr = std::get<std::shared_ptr<NaryOp<T, T, T>>>(_ops.back()).get();
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            else if (llvm::isa<mlir::ComplexType>(eleLhsTy) &&
                llvm::isa<mlir::Float64Type>(eleLhsTy.cast<mlir::ComplexType>().getElementType()) &&
                llvm::isa<mlir::ComplexType>(eleRhsTy) &&
                llvm::isa<mlir::Float64Type>(eleRhsTy.cast<mlir::ComplexType>().getElementType()))
            {
                using ptTy = std::complex<double>;
                using T = View::View<ptTy, View::DCCSC3D>;
                using op_t = Op<SubFunctor<ptTy>, T, T, T>;
                _ops.push_back(std::make_unique<op_t>(SubFunctor<ptTy>()));
                auto* ptr = std::get<std::shared_ptr<NaryOp<T, T, T>>>(_ops.back()).get();
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            else
            {
                throw std::logic_error("sub op not implemented.");
            }
        }
        #ifdef QUICC_HAS_CUDA_BACKEND
        else
        {
            using namespace QuICC::Pointwise::Cuda;
            using namespace QuICC::Pointwise;
            if (llvm::isa<mlir::Float64Type>(eleLhsTy) &&
                llvm::isa<mlir::Float64Type>(eleRhsTy))
            {
                using ptTy = double;
                using T = R_DCCSC3D_t;
                using op_t = Op<SubFunctor<ptTy>, T, T, T>;
                _ops.push_back(std::make_unique<op_t>(SubFunctor<ptTy>()));
                auto* ptr = std::get<std::shared_ptr<NaryOp<T, T, T>>>(_ops.back()).get();
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            else if (llvm::isa<mlir::ComplexType>(eleLhsTy) &&
                llvm::isa<mlir::Float64Type>(eleLhsTy.cast<mlir::ComplexType>().getElementType()) &&
                llvm::isa<mlir::ComplexType>(eleRhsTy) &&
                llvm::isa<mlir::Float64Type>(eleRhsTy.cast<mlir::ComplexType>().getElementType()))
            {
                using ptTy = std::complex<double>;
                using T = View::View<ptTy, View::DCCSC3D>;
                using op_t = Op<SubFunctor<ptTy>, T, T, T>;
                _ops.push_back(std::make_unique<op_t>(SubFunctor<ptTy>()));
                auto* ptr = std::get<std::shared_ptr<NaryOp<T, T, T>>>(_ops.back()).get();
                assert(ptr != nullptr);
                _thisArr[index] = ptr;
            }
            else
            {
                throw std::logic_error("sub op not implemented.");
            }

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
