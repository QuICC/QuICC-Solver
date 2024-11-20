#include "Graph/OpsMap.hpp"

#include "Graph/Shims/MlirShims.hpp"

namespace QuICC
{
namespace Graph
{

MapOps::MapOps(mlir::ModuleOp module,
      const PhysicalParameters<double> physParams,
      const std::shared_ptr<Memory::memory_resource> mem) :
        _mem(mem), _physParams(physParams)
         {
    #ifdef QUICC_HAS_CUDA_BACKEND
    // Check memory space
    {
      QuICC::Memory::MemBlock<std::size_t> test(1, _mem.get());
      _isCpu = !QuICC::Cuda::isDeviceMemory(test.data());
    }
    #endif

    #ifdef QUICC_MPI
    _commFTAL = std::make_shared<Transpose::Mpi::Comm<std::complex<double>>>();
    _commALFT = std::make_shared<Transpose::Mpi::Comm<std::complex<double>>>();
    _commALJW = std::make_shared<Transpose::Mpi::Comm<std::complex<double>>>();
    _commJWAL = std::make_shared<Transpose::Mpi::Comm<std::complex<double>>>();
    #endif

    using namespace mlir;
    Dialect *quiccirDialect = module->getContext()->getLoadedDialect("quiccir");
    mlir::WalkResult result = module->walk([&](mlir::Operation* op) {
#ifndef NDEBUG
      if (op->getDialect() == quiccirDialect) {
        llvm::outs() << "visiting: " << op->getName() << '\n';
      }
#endif
      if (auto frPrj = dyn_cast<mlir::quiccir::FrPOp>(op)) {
        setFourierPrj(frPrj);
      }
      else if (auto frInt = dyn_cast<mlir::quiccir::FrIOp>(op)) {
        setFourierInt(frInt);
      }
      else if (auto alPrj = dyn_cast<mlir::quiccir::AlPOp>(op)) {
        setALegendrePrj(alPrj);
      }
      else if (auto alInt = dyn_cast<mlir::quiccir::AlIOp>(op)) {
        setALegendreInt(alInt);
      }
      else if (auto jwPrj = dyn_cast<mlir::quiccir::JWPOp>(op)) {
        setWorlandPrj(jwPrj);
      }
      else if (auto jwInt = dyn_cast<mlir::quiccir::JWIOp>(op)) {
        setWorlandInt(jwInt);
      }
      else if (auto add = dyn_cast<mlir::quiccir::AddOp>(op)) {
        setAdd(add);
      }
      else if (auto sub = dyn_cast<mlir::quiccir::SubOp>(op)) {
        setSub(sub);
      }
      else if (auto cross = dyn_cast<mlir::quiccir::CrossOp>(op)) {
        setCross(cross);
      }
      else if (auto dot = dyn_cast<mlir::quiccir::DotOp>(op)) {
        setDot(dot);
      }
      else if (auto mulConst = dyn_cast<mlir::quiccir::MulConstOp>(op)) {
        setMulConst(mulConst);
      }
      else if (auto tran = dyn_cast<mlir::quiccir::TransposeOp>(op)) {
        setTranspose(tran);
      }
      // return deallocateBuffers(op);
    //   if (failed(deallocateBuffers(op)))
    //     return WalkResult::interrupt();
      return WalkResult::advance();
    });

}


} // namespace Graph
} // namespace QuICC

