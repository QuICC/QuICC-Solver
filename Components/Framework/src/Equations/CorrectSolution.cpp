/**
 * @file CorrectSolution.cpp
 * @brief CorrectSolution function
 */

// System includes
//

// Project includes
//
#include "QuICC/Equations/CorrectSolution.hpp"

namespace QuICC {

namespace Equations {

   std::vector<std::tuple<MHDComplex,int,int>> correctSolution(const IFieldEquation& eq, FieldComponents::Spectral::Id compId, const std::vector<std::tuple<MHDVariant,int,int,int>>& corr, const int matIdx, const int start)
   {
      // matIdx is the index of the slowest varying direction with a single RHS
      if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         details::CorrectSolutionFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS> func(eq, corr, matIdx);
         return func.corr;
      }
      // matIdx is the index of the slowest varying direction with multiple RHS
      else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
      {
         details::CorrectSolutionFunctor<CouplingIndexType::SLOWEST_MULTI_RHS> func(eq, corr, matIdx);
         return func.corr;
      }
      // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
      else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::MODE)
      {
         details::CorrectSolutionFunctor<CouplingIndexType::MODE> func(eq, corr, matIdx);
         return func.corr;
      }
      // There is a single matrix
      else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SINGLE)
      {
         details::CorrectSolutionFunctor<CouplingIndexType::SINGLE> func(eq, corr, matIdx);
         return func.corr;
      }
      else
      {
         throw std::logic_error("Unknown coupling");
      }
   }

} // Equations
} // QuICC
