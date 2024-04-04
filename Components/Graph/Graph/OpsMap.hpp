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

// al ops
#include "QuICC/Polynomial/ALegendre/Plm.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
#include "ViewOps/ALegendre/Builder.hpp"

#include "ViewOps/Pointwise/Cpu/Pointwise.hpp"

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
    std::shared_ptr<UnaryOp<R_DCCSC3D_t, R_DCCSC3D_t>>,
    std::shared_ptr<UnaryOp<C_S1CLCSC3D_t, C_DCCSC3D_t>>
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
      const std::shared_ptr<Memory::memory_resource> mem,
      const std::array<std::uint32_t, 3> physDims,
      const std::array<std::uint32_t, 3> modsDims);

    /// @brief return void pointers to ops
    std::vector<void*> getThisArr(){return _thisArr;}

};


MapOps::MapOps(mlir::ModuleOp module,
      const std::shared_ptr<Memory::memory_resource> mem,
      const std::array<std::uint32_t, 3> physDims,
      const std::array<std::uint32_t, 3> modsDims) : _mem(mem) {
    using namespace mlir;
    Dialect *quiccirDialect = module->getContext()->getLoadedDialect("quiccir");
    mlir::WalkResult result = module->walk([&](mlir::Operation* op) {
#ifndef NDEBUG
      if (op->getDialect() == quiccirDialect) {
        llvm::outs() << "visiting: " << op->getName() << '\n';
      }
#endif
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
        auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
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
        auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
        assert(ptr != nullptr);
        _thisArr[index] = ptr;

      }
      else if (auto alInt = dyn_cast<mlir::quiccir::AlIOp>(op)) {
        using namespace QuICC::Transform::Quadrature;
        using Tin = C_DCCSC3D_t;
        using Tout = C_S1CLCSC3D_t;
        using Top = QuICC::View::View<double, QuICC::View::S1CLCSC3DJIK>;
        using backend_t = Cpu::ImplOp<Tout, Tin, Top>;
        using op_t = Op<Tout, Tin, Top, backend_t>;
        _ops.push_back(std::make_unique<op_t>(mem));
        // Get index from MLIR source
        std::uint64_t index = alInt.getImplptr().value();
        if (index >= _thisArr.size()) {
          _thisArr.resize(index+1);
        }
        auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
        // Setup operator
        auto alIntOp = dynamic_cast<op_t*>(ptr);
        constexpr size_t rank = 3;
        /// dim 0 - L  - harmonic degree
        /// dim 1 - Nl - longitudinal points
        /// dim 2 - M  - harmonic order
        std::array<std::uint32_t, rank> dims {modsDims[1], physDims[1], physDims[2]};
        std::vector<std::uint32_t> layers;
        // Dense operator \todo generalize for distributed op
        for (std::size_t i = 0; i < dims[2]; ++i) {
          layers.push_back(i);
        }
        alIntOp->allocOp(dims, layers);
        // Set grid \todo set once per operator kind
        Internal::Array igrid;
        Internal::Array iweights;
        ::QuICC::Polynomial::Quadrature::LegendreRule quad;
        quad.computeQuadrature(igrid, iweights, physDims[1]);
        // Populate op
        auto opView = alIntOp->getOp();
        using namespace QuICC::Transform::ALegendre;
        builder<Top, Polynomial::ALegendre::Plm, Internal::Array::Scalar, 0>(opView, igrid, iweights);
        // Add to thisArr
        assert(ptr != nullptr);
        _thisArr[index] = ptr;

      }
      else if (auto add = dyn_cast<mlir::quiccir::AddOp>(op)) {
        using namespace QuICC::Pointwise::Cpu;
        using namespace QuICC::Pointwise;
        using T = C_DCCSC3D_t;
        using op_t = Op<AddFunctor<std::complex<double>>, T, T, T>;
        _ops.push_back(std::make_unique<op_t>(AddFunctor<std::complex<double>>()));
        // get index from MLIR source
        std::uint64_t index = add.getImplptr().value();
        if (index >= _thisArr.size()) {
          _thisArr.resize(index+1);
        }
        auto* ptr = std::get<std::shared_ptr<NaryOp<C_DCCSC3D_t, C_DCCSC3D_t, C_DCCSC3D_t>>>(_ops.back()).get();
        assert(ptr != nullptr);
        _thisArr[index] = ptr;
      }
      else if (auto tran = dyn_cast<mlir::quiccir::TransposeOp>(op)) {
        using namespace QuICC::Transpose::Cpu;
        using namespace QuICC::Transpose;
        using Tin = C_DCCSC3D_t;
        using Tout = C_DCCSC3D_t;
        using op_t = Op<Tout, Tin, p021_t>;
        _ops.push_back(std::make_unique<op_t>());
        // get index from MLIR source
        std::uint64_t index = tran.getImplptr().value();
        if (index >= _thisArr.size()) {
          _thisArr.resize(index+1);
        }
        auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
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

