/**
 * @file CorrectSolution.hpp
 * @brief CorrectSolution function
 */

#ifndef QUICC_EQUATIONS_CORRECTSOLUTION_HPP
#define QUICC_EQUATIONS_CORRECTSOLUTION_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/IEquation.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/Equations/details/CorrectSolutionFunctor.hpp"
#include "QuICC/Equations/details/CorrectSolutionFunctorM.hpp"
#include "QuICC/Equations/details/CorrectSolutionFunctorS.hpp"
#include "QuICC/Equations/details/CorrectSolutionFunctorSMR.hpp"
#include "QuICC/Equations/details/CorrectSolutionFunctorSSR.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "Arithmetics/Basic.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Correct solution based on constraint correction
    *
    * @param eq         Equation to work on
    * @param compId     Component ID
    * @param storage    Storage for the equation values
    * @param matIdx     Index of the given data
    * @param start      Start index for the storage
    */
   template <typename TData> void correctSolution(const IFieldEquation& eq, FieldComponents::Spectral::Id compId, const std::vector<std::tuple<MHDVariant,int,int,int>>& corr, TData& storage, const int matIdx, const int start);

   /**
    * @brief Remap constraint corrections
    *
    * @param eq         Equation to work on
    * @param compId     Component ID
    * @param matIdx     Index of the given data
    * @param start      Start index for the storage
    */
   std::vector<std::tuple<MHDComplex,int,int>> correctSolution(const IFieldEquation& eq, FieldComponents::Spectral::Id compId, const std::vector<std::tuple<MHDVariant,int,int,int>>& corr, const int matIdx, const int start);

   template <typename TData> void correctSolution(const IFieldEquation& eq, FieldComponents::Spectral::Id compId, const std::vector<std::tuple<MHDVariant,int,int,int>>& corr, TData& storage, const int matIdx, const int start)
   {
      // matIdx is the index of the slowest varying direction with a single RHS
      if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         details::CorrectSolutionFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS> func(eq, corr, matIdx);
         func.apply(storage, start);
      }
      // matIdx is the index of the slowest varying direction with multiple RHS
      else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
      {
         details::CorrectSolutionFunctor<CouplingIndexType::SLOWEST_MULTI_RHS> func(eq, corr, matIdx);
         func.apply(storage, start);
      }
      // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
      else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::MODE)
      {
         details::CorrectSolutionFunctor<CouplingIndexType::MODE> func(eq, corr, matIdx);
         func.apply(storage, start);
      }
      // There is a single matrix
      else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SINGLE)
      {
         details::CorrectSolutionFunctor<CouplingIndexType::SINGLE> func(eq, corr, matIdx);
         func.apply(storage, start);
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_CORRECTSOLUTION_HPP
