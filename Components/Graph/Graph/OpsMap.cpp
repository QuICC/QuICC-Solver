#include "Graph/OpsMap.hpp"

#include "Graph/BackendsMap.hpp"
#include "Graph/Shims/MlirShims.hpp"

namespace QuICC
{
namespace Graph
{

MapOps::MapOps(mlir::ModuleOp module,
      const std::shared_ptr<Memory::memory_resource> mem,
      const std::array<std::uint32_t, 3> physDims,
      const std::array<std::uint32_t, 3> modsDims) :
      _mem(mem), _physDims(physDims), _modsDims(modsDims) {
    #ifdef QUICC_HAS_CUDA_BACKEND
    // Check memory space
    {
      QuICC::Memory::MemBlock<std::size_t> test(1, _mem.get());
      _isCpu = !QuICC::Cuda::isDeviceMemory(test.data());
    }
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
      else if (auto alInt = dyn_cast<mlir::quiccir::AlIOp>(op)) {
        setALegendreInt(alInt);
      }
      else if (auto jwInt = dyn_cast<mlir::quiccir::JWIOp>(op)) {
        setWorlandInt(jwInt);
      }
      // return deallocateBuffers(op);
    //   if (failed(deallocateBuffers(op)))
    //     return WalkResult::interrupt();
      return WalkResult::advance();
    });

}


} // namespace Graph
} // namespace QuICC

