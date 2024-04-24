#pragma once

#include <vector>
#include <cassert>
#include <iostream>
#include <Quiccir/IR/QuiccirDialect.h>
#include <Quiccir/IR/QuiccirOps.h>
#include <Quiccir/Pipelines/Passes.h>
#include <mlir/IR/BuiltinDialect.h>
#include <mlir/IR/Operation.h>

#include "Graph/BackendsMap.hpp"
#include "Graph/Shims/MlirShims.hpp"
#include "Graph/Types.hpp"

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

/// @brief map mlir ops to QuICC operators
class MapOps
{
  private:
    /// @brief store pointers for passing into MLIR
    std::vector<void*> _thisArr;
    /// @brief store for RAII
    std::vector<varOp_t> _ops;
    /// @brief memory resource for internal op allocation
    std::shared_ptr<Memory::memory_resource> _mem;
    /// @brief is the memory space cpu or gpu
    /// to be replaced by a per-operator-attribute
    bool _isCpu = true;
  public:
    /// @brief empty constructor
    MapOps() = default;
    /// @brief constructor from mlir module
    /// @param module mlir module
    /// @param mem memory resource to pass to operators
    /// @param physDims physical dimensions of the problem
    /// @param modsDims spectral dimensions of the problem
    MapOps(mlir::ModuleOp module,
      const std::shared_ptr<Memory::memory_resource> mem,
      const std::array<std::uint32_t, 3> physDims,
      const std::array<std::uint32_t, 3> modsDims);

    /// @brief return void pointers to ops
    std::vector<void*> getThisArr(){return _thisArr;}

  private:
    /// @brief map Fourier projectors
    /// @param op
    void setFourierPrj(mlir::quiccir::FrPOp op);
    /// @brief map Fourier integrators
    /// @param op
    void setFourierInt(mlir::quiccir::FrIOp op);

};

} // namespace Graph
} // namespace QuICC

