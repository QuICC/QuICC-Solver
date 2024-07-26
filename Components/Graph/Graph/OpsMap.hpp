/**
 * @file OPsMap.hpp
 * @brief map mlir operators to QuICC
 */
#pragma once

// External includes
//
#include <vector>
#include <Quiccir/IR/QuiccirDialect.h>
#include <Quiccir/IR/QuiccirOps.h>
#include <Quiccir/Pipelines/Passes.h>
#include <mlir/IR/BuiltinDialect.h>
#include <mlir/IR/Operation.h>

// Project includes
//
#include "Graph/Types.hpp"
#include "Memory/Memory.hpp"
#ifdef QUICC_MPI
#include "ViewOps/Transpose/Mpi/Comm.hpp"
#endif

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
    /// comm
    #ifdef QUICC_MPI
    std::shared_ptr<Transpose::Mpi::Comm<std::complex<double>>> _commFTAL;
    std::shared_ptr<Transpose::Mpi::Comm<std::complex<double>>> _commALFT;
    std::shared_ptr<Transpose::Mpi::Comm<std::complex<double>>> _commALJW;
    std::shared_ptr<Transpose::Mpi::Comm<std::complex<double>>> _commJWAL;
    #endif
  public:
    /// @brief empty constructor
    MapOps() = default;
    /// @brief constructor from mlir module
    /// @param module mlir module
    /// @param mem memory resource to pass to operators
    MapOps(mlir::ModuleOp module,
      const std::shared_ptr<Memory::memory_resource> mem);

    /// @brief return void pointers to ops
    std::vector<void*> getThisArr(){return _thisArr;}

  private:
    /// @brief map Fourier projectors
    /// @param op
    void setFourierPrj(mlir::quiccir::FrPOp op);
    /// @brief map Fourier integrators
    /// @param op
    void setFourierInt(mlir::quiccir::FrIOp op);
    /// @brief map ALegendre projectors
    /// @param op
    void setALegendrePrj(mlir::quiccir::AlPOp op);
    /// @brief map ALegendre integrators
    /// @param op
    void setALegendreInt(mlir::quiccir::AlIOp op);
    /// @brief map ALegendre projectors
    /// @param op
    void setWorlandPrj(mlir::quiccir::JWPOp op);
    /// @brief map Worland integrators
    /// @param op
    void setWorlandInt(mlir::quiccir::JWIOp op);
    /// @brief map pointwise addition
    /// @param op
    void setAdd(mlir::quiccir::AddOp op);
    /// @brief map pointwise subtraction
    /// @param op
    void setSub(mlir::quiccir::SubOp op);
    /// @brief map transpose
    /// @param op
    void setTranspose(mlir::quiccir::TransposeOp op);
};

} // namespace Graph
} // namespace QuICC

