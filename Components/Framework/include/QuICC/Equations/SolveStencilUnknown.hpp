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
#include "QuICC/Equations/Dispatchers.hpp"
#include "QuICC/Solver/SparseSolver.hpp"
#include "QuICC/SparseSolvers/SparseLinearSolverTools.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Solve for galerkin unknown using the stencil
    *
    * @param compId     Component ID
    * @param storage    Storage for the equation values
    * @param matIdx     Index of the given data
    * @param start      Start index for the storage
    */
   template <typename TField, typename TData> void solveStencilUnknown(const Resolution& res, const CouplingInformation& cinfo, std::size_t fieldName, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const Model::IModelBackend& backend, const std::map<std::size_t, std::size_t>& bcIds, const std::map<std::size_t, NonDimensional::SharedINumber>& eqParams);

   template <typename TField, typename TData> void solveStencilUnknown(const Resolution& res, const CouplingInformation& cinfo, std::size_t fieldName, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const Model::IModelBackend& backend, const std::map<std::size_t, std::size_t>& bcIds, const std::map<std::size_t, NonDimensional::SharedINumber>& eqParams)
   {
      using TmpDataType = typename std::conditional<std::is_same_v<TData, DecoupledZMatrix>, DecoupledZMatrix, Eigen::Matrix<typename Arithmetics::GetScalarType<TData>::ScalarType, Eigen::Dynamic, Eigen::Dynamic>>::type;

      // Create temporary storage for tau data
      TmpDataType tmp(cinfo.tauN(matIdx), cinfo.rhsCols(matIdx));
      Equations::copyUnknown(res, cinfo, field, compId, tmp, matIdx, 0, false, true, true);
      TmpDataType rhs(cinfo.galerkinN(matIdx), cinfo.rhsCols(matIdx));
      if(res.sim().ss().has(SpatialScheme::Feature::SpectralMatrix2D))
      {
         Arithmetics::setTopBlock(rhs, 0, cinfo.galerkinN(matIdx), res.sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL), cinfo.galerkinShift(matIdx, 0), tmp);
      }
      else
      {
         Arithmetics::setTopBlock(rhs, 0, cinfo.galerkinN(matIdx), tmp);
      }

      // Get a restricted stencil matrix
      SparseMatrix stencil(cinfo.galerkinN(matIdx),cinfo.galerkinN(matIdx));
      dispatchGalerkinStencil(fieldName, compId, stencil, matIdx, res, true, backend, cinfo, bcIds, eqParams);
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
      TmpDataType lhs(cinfo.galerkinN(matIdx), cinfo.rhsCols(matIdx));
      Solver::details::solveWrapper(lhs, solver, rhs);
      Arithmetics::setTopBlock(storage, start, cinfo.galerkinN(matIdx), lhs);
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_SOLVESTENCILUNKNOWN_HPP
