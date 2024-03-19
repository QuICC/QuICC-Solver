#pragma once

#include <vector>
#include <cassert>
#include <iostream>
#include <Quiccir/IR/QuiccirDialect.h>
#include <Quiccir/IR/QuiccirOps.h>
#include <Quiccir/Pipelines/Passes.h>
#include <mlir/IR/BuiltinDialect.h>
#include <mlir/IR/Operation.h>

#include "Operator/Nary.hpp"
#include "Operator/Unary.hpp"
#include "Graph/MlirShims.hpp"

namespace QuICC
{
namespace Graph
{


using QuICC::Operator::NaryOp;
using QuICC::Operator::UnaryOp;

/// @brief encode unary ops
using varOp_t = std::variant<
    // std::shared_ptr<NaryOp<R_VB_t, R_VB_t>>,
    std::shared_ptr<NaryOp<C_DCCSC3D_t, C_DCCSC3D_t>>,
    std::shared_ptr<NaryOp<R_DCCSC3D_t, R_DCCSC3D_t>>,
    // std::shared_ptr<NaryOp<R_VB_t, R_VB_t, R_VB_t>>,
    std::shared_ptr<NaryOp<C_DCCSC3D_t, C_DCCSC3D_t, C_DCCSC3D_t>>,
    std::shared_ptr<NaryOp<R_DCCSC3D_t, R_DCCSC3D_t, R_DCCSC3D_t>>,
    std::shared_ptr<UnaryOp<C_DCCSC3D_t, C_DCCSC3D_t>>,
    std::shared_ptr<UnaryOp<R_DCCSC3D_t, C_DCCSC3D_t>>,
    std::shared_ptr<UnaryOp<C_DCCSC3D_t, R_DCCSC3D_t>>,
    std::shared_ptr<UnaryOp<R_DCCSC3D_t, R_DCCSC3D_t>>
>;


class MapOps
{
  private:
    /// @brief store pointers for passing into MLIR
    std::vector<void*> _thisArr;
    /// @brief store for RAII
    std::vector<varOp_t> _ops;
    /// @brief memory resource for internal op allocation
    const std::shared_ptr<Memory::memory_resource> _mem;
  public:
    MapOps(mlir::ModuleOp module,
      const std::shared_ptr<Memory::memory_resource> mem);

    /// @brief return void pointers to ops
    std::vector<void*> getThisArr(){return _thisArr;}

};


MapOps::MapOps(mlir::ModuleOp module,
      const std::shared_ptr<Memory::memory_resource> mem) : _mem(mem) {
    using namespace mlir;
    Dialect *quiccirDialect = module->getContext()->getLoadedDialect("quiccir");
    mlir::WalkResult result = module->walk([&](mlir::Operation* op) {
      if (op->getDialect() == quiccirDialect) {
        llvm::outs() << op->getName() << '\n';
      }
      if (auto frPrj = dyn_cast<mlir::quiccir::FrPOp>(op)) {
        using namespace QuICC::Transform::Fourier;
        using backend_t = QuICC::Graph::viewCpu_t;
        using Tin = C_DCCSC3D_t;
        using Tout = R_DCCSC3D_t;
        using backendFft_t = Fft_t<backend_t, Tout, Tin>;
        using backendDiff_t = MixedDiff_t<backend_t, Tin, 0, bwd_t,
          QuICC::Transform::Fourier::none_m>;
        using op_t = Mixed::Projector::DOp<Tout, Tin, backendFft_t,
        backendDiff_t>;
        _ops.push_back(std::make_unique<op_t>(mem));
        // get index from MLIR source
        std::uint64_t index = frPrj.getImplptr().value();
        if (index >= _thisArr.size()) {
          _thisArr.resize(index+1);
        }
        auto* ptr = std::get<std::shared_ptr<UnaryOp<R_DCCSC3D_t, C_DCCSC3D_t>>>(_ops.back()).get();
        assert(ptr != nullptr);
        _thisArr[index] = ptr;

      }
      else if (auto frInt = dyn_cast<mlir::quiccir::FrIOp>(op)) {
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
        // get index from MLIR source
        std::uint64_t index = frInt.getImplptr().value();
        if (index >= _thisArr.size()) {
          _thisArr.resize(index+1);
        }
        auto* ptr = std::get<std::shared_ptr<UnaryOp<C_DCCSC3D_t, R_DCCSC3D_t>>>(_ops.back()).get();
        assert(ptr != nullptr);
        _thisArr[index] = ptr;

      }
      // return deallocateBuffers(op);
    //   if (failed(deallocateBuffers(op)))
    //     return WalkResult::interrupt();
      return WalkResult::advance();
    });

}


} // namespace Graph
} // namespace QuICC

