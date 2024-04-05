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

// jw ops
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandRule.hpp"
#include "ViewOps/Worland/Tags.hpp"
#include "ViewOps/Worland/Builder.hpp"
#include "DenseSM/Worland/Operator.hpp"

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
        // Get index from MLIR source
        std::uint64_t index = frPrj.getImplptr().value();
        if (index >= _thisArr.size()) {
          _thisArr.resize(index+1, nullptr);
        }
        if (_thisArr[index] == nullptr) {
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
          auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
          assert(ptr != nullptr);
          _thisArr[index] = ptr;
        }
        else {
          #ifndef NDEBUG
          std::cout << "operator already allocated\n";
          #endif
        }
      }
      else if (auto frInt = dyn_cast<mlir::quiccir::FrIOp>(op)) {
        // Get index from MLIR source
        std::uint64_t index = frInt.getImplptr().value();
        if (index >= _thisArr.size()) {
          _thisArr.resize(index+1, nullptr);
        }
        if (_thisArr[index] == nullptr) {
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
          auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
          assert(ptr != nullptr);
          _thisArr[index] = ptr;
        }
        else {
          #ifndef NDEBUG
          std::cout << "operator already allocated\n";
          #endif
        }
      }
      else if (auto alInt = dyn_cast<mlir::quiccir::AlIOp>(op)) {
        // Get index from MLIR source
        std::uint64_t index = alInt.getImplptr().value();
        if (index >= _thisArr.size()) {
          _thisArr.resize(index+1, nullptr);
        }
        if (_thisArr[index] == nullptr) {
          using namespace QuICC::Transform::Quadrature;
          using Tin = C_DCCSC3D_t;
          using Tout = C_S1CLCSC3D_t;
          using Top = QuICC::View::View<double, QuICC::View::S1CLCSC3DJIK>;
          using backend_t = Cpu::ImplOp<Tout, Tin, Top>;
          using op_t = Op<Tout, Tin, Top, backend_t>;
          _ops.push_back(std::make_unique<op_t>(mem));
          auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
          // Setup operator
          auto alIntOp = dynamic_cast<op_t*>(ptr);
          constexpr size_t rank = 3;
          /// dim 0 - L  - harmonic degree
          /// dim 1 - Ntheta - colatitudinal points
          /// dim 2 - M  - harmonic order
          std::array<std::uint32_t, rank> dims {modsDims[1], physDims[1], modsDims[0]};
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
          quad.computeQuadrature(igrid, iweights, dims[1]);
          // Populate op
          auto opView = alIntOp->getOp();
          using namespace QuICC::Transform::ALegendre;
          builder<Top, Polynomial::ALegendre::Plm, Internal::Array::Scalar, 0>(opView, igrid, iweights);
          // Add to thisArr
          assert(ptr != nullptr);
          _thisArr[index] = ptr;
        }
        else {
          #ifndef NDEBUG
          std::cout << "operator already allocated\n";
          #endif
        }
      }
      else if (auto jwInt = dyn_cast<mlir::quiccir::JWIOp>(op)) {
        // Get index from MLIR source
        std::uint64_t index = jwInt.getImplptr().value();
        if (index >= _thisArr.size()) {
          _thisArr.resize(index+1, nullptr);
        }
        if (_thisArr[index] == nullptr) {
          using namespace QuICC::Transform::Quadrature;
          using Tin = C_DCCSC3D_t;
          using Tout = C_DCCSC3D_t;
          using Top = QuICC::View::View<double, QuICC::View::CSL3DJIK>;
          using backend_t = Cpu::ImplOp<Tout, Tin, Top>;
          using op_t = Op<Tout, Tin, Top, backend_t>;
          _ops.push_back(std::make_unique<op_t>(mem));
          auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
          // Setup operator
          auto jwIntOp = dynamic_cast<op_t*>(ptr);
          constexpr size_t rank = 3;
          // dim 0 - N  - radial modes
          // dim 1 - Nr - radial points
          // dim 2 - L  - harmonic degree
          std::array<std::uint32_t, rank> dims {modsDims[2], physDims[2], modsDims[1]};
          std::vector<std::uint32_t> layers;
          /// Dense operator \todo generalize for distributed op
          for (std::size_t i = 0; i < dims[2]; ++i) {
            layers.push_back(i);
          }
          jwIntOp->allocOp(dims, layers);
          // Set grid \todo set once per operator kind
          Internal::Array igrid;
          Internal::Array iweights;
          ::QuICC::Polynomial::Quadrature::WorlandRule quad;
          quad.computeQuadrature(igrid, iweights, dims[1]);
          // Populate op
          auto opView = jwIntOp->getOp();
          using namespace QuICC::Polynomial::Worland;
          using namespace QuICC::Transform::Worland;
          Builder<Top, QuICC::DenseSM::Worland::Operator<Wnl>, fwd_t> tBuilder;
          tBuilder.compute(opView, igrid, iweights);
          // Add to thisArr
          assert(ptr != nullptr);
          _thisArr[index] = ptr;
        }
        else {
          #ifndef NDEBUG
          std::cout << "operator already allocated\n";
          #endif
        }
      }
      else if (auto add = dyn_cast<mlir::quiccir::AddOp>(op)) {
        // Get index from MLIR source
        std::uint64_t index = add.getImplptr().value();
        if (index >= _thisArr.size()) {
          _thisArr.resize(index+1, nullptr);
        }
        if (_thisArr[index] == nullptr) {
          using namespace QuICC::Pointwise::Cpu;
          using namespace QuICC::Pointwise;
          using T = C_DCCSC3D_t;
          using op_t = Op<AddFunctor<std::complex<double>>, T, T, T>;
          _ops.push_back(std::make_unique<op_t>(AddFunctor<std::complex<double>>()));
          auto* ptr = std::get<std::shared_ptr<NaryOp<C_DCCSC3D_t, C_DCCSC3D_t, C_DCCSC3D_t>>>(_ops.back()).get();
          assert(ptr != nullptr);
          _thisArr[index] = ptr;
        }
        else {
          #ifndef NDEBUG
          std::cout << "operator already allocated\n";
          #endif
        }
      }
      else if (auto sub = dyn_cast<mlir::quiccir::SubOp>(op)) {
        // Get index from MLIR source
        std::uint64_t index = sub.getImplptr().value();
        if (index >= _thisArr.size()) {
          _thisArr.resize(index+1, nullptr);
        }
        if (_thisArr[index] == nullptr) {
          using namespace QuICC::Pointwise::Cpu;
          using namespace QuICC::Pointwise;
          using T = C_DCCSC3D_t;
          using op_t = Op<SubFunctor<std::complex<double>>, T, T, T>;
          _ops.push_back(std::make_unique<op_t>(SubFunctor<std::complex<double>>()));
          auto* ptr = std::get<std::shared_ptr<NaryOp<C_DCCSC3D_t, C_DCCSC3D_t, C_DCCSC3D_t>>>(_ops.back()).get();
          assert(ptr != nullptr);
          _thisArr[index] = ptr;
        }
        else {
          #ifndef NDEBUG
          std::cout << "operator already allocated\n";
          #endif
        }
      }
      else if (auto tran = dyn_cast<mlir::quiccir::TransposeOp>(op)) {
        // Get index from MLIR source
        std::uint64_t index = tran.getImplptr().value();
        if (index >= _thisArr.size()) {
          _thisArr.resize(index+1, nullptr);
        }
        if (_thisArr[index] == nullptr) {
          using namespace QuICC::Transpose::Cpu;
          using namespace QuICC::Transpose;
          // check type attribute
          using Tin = C_DCCSC3D_t;
          using Tout = C_DCCSC3D_t;
          // check perm attribute
          auto perm = tran.getPermutation();
          assert(perm[0] == 2);
          assert(perm[1] == 0);
          assert(perm[2] == 1);
          using op_t = Op<Tout, Tin, p201_t>;
          _ops.push_back(std::make_unique<op_t>());
          auto* ptr = std::get<std::shared_ptr<UnaryOp<Tout, Tin>>>(_ops.back()).get();
          assert(ptr != nullptr);
          _thisArr[index] = ptr;
        }
        else {
          #ifndef NDEBUG
          std::cout << "operator already allocated\n";
          #endif
        }
      }
      // return deallocateBuffers(op);
    //   if (failed(deallocateBuffers(op)))
    //     return WalkResult::interrupt();
      return WalkResult::advance();
    });

}


} // namespace Graph
} // namespace QuICC

