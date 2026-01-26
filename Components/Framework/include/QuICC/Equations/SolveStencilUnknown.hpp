/**
 *
 * @file SolveStencilUnknown.hpp
 * @brief Base for the implementation of a scalar equation
 */

#ifndef QUICC_EQUATIONS_SOLVESTENCILUNKNOWN_HPP
#define QUICC_EQUATIONS_SOLVESTENCILUNKNOWN_HPP

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Arithmetics/Utility.hpp"
#include "Eigen/src/Core/util/Constants.h"
#include "Types/Typedefs.hpp"
#include "QuICC/Equations/CopyUnknown.hpp"
#include "QuICC/Solver/SparseSolver.hpp"
#include "QuICC/SparseSolvers/SparseLinearSolverTools.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Solve for galerkin unknown using the stencil
    *
    * @param eq         Equation to work on
    * @param compId     Component ID
    * @param storage    Storage for the equation values
    * @param matIdx     Index of the given data
    * @param start      Start index for the storage
    */
   template <typename TEquation, typename TData> void solveStencilUnknown(const TEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   template <typename TEquation, typename TData> void solveStencilUnknown(const TEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      using TmpDataType = typename std::conditional<std::is_same_v<TData, DecoupledZMatrix>, DecoupledZMatrix, Eigen::Matrix<typename Arithmetics::GetScalarType<TData>::ScalarType, Eigen::Dynamic, Eigen::Dynamic>>::type;

      const auto& info = eq.couplingInfo(compId);

      // Create temporary storage for tau data
      TmpDataType tmp(info.tauN(matIdx), info.rhsCols(matIdx));
      std::visit(
            [&](auto&& p)
            {
               Equations::copyUnknown(eq, p->dom(0).perturbation(), compId, tmp, matIdx, 0, false, true);
            }, eq.spUnknown());
      TmpDataType rhs(info.galerkinN(matIdx), info.rhsCols(matIdx));
      if(eq.res().sim().ss().has(SpatialScheme::Feature::SpectralMatrix2D))
      {
         Arithmetics::setTopBlock(rhs, 0, info.galerkinN(matIdx), eq.res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL), info.galerkinShift(matIdx, 0), tmp);
      }
      else
      {
         Arithmetics::setTopBlock(rhs, 0, info.galerkinN(matIdx), tmp);
      }

      // Get a restricted stencil matrix
      SparseMatrix stencil(info.galerkinN(matIdx),info.galerkinN(matIdx));
      eq.dispatchGalerkinStencil(compId, stencil, matIdx, eq.res(), info.couplingTools().getIndexes(eq.res(), matIdx), true);
      stencil.makeCompressed();

      // Check that square stencil was generated. Setup is wrong if matrix is not square
      if(stencil.rows() != stencil.cols())
      {
      	throw std::logic_error("Stencil setup is wrong and did not produce a square matrix");
      }

      // Create solver and factorize stencil
      Framework::Selector::SparseSolver<SparseMatrix> solver;
      solver.compute(stencil);
      // Safety assert for successful factorisation
      if(solver.info() != Eigen::Success)
      {
         throw std::logic_error("Stencil factorization for initial solution failed!");
      }

      // solve for galerkin expansion
      TmpDataType lhs(info.galerkinN(matIdx), info.rhsCols(matIdx));
      Solver::details::solveWrapper(lhs, solver, rhs);
      Arithmetics::setTopBlock(storage, start, info.galerkinN(matIdx), lhs);
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_SOLVESTENCILUNKNOWN_HPP
