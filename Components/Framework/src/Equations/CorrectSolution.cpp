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

   std::vector<std::tuple<MHDComplex,int,int>> correctSolution(const Resolution& res, const CouplingInformation& cinfo, FieldComponents::Spectral::Id compId, const std::vector<std::tuple<MHDVariant,int,int,int>>& corr, const int matIdx, const int start)
   {
      // matIdx is the index of the slowest varying direction with a single RHS
      if(cinfo.indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         details::CorrectSolutionFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS> func(res, cinfo, corr, matIdx);
         return func.corr;
      }
      // matIdx is the index of the slowest varying direction with multiple RHS
      else if(cinfo.indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
      {
         details::CorrectSolutionFunctor<CouplingIndexType::SLOWEST_MULTI_RHS> func(res, cinfo, corr, matIdx);
         return func.corr;
      }
      // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
      else if(cinfo.indexType() == CouplingIndexType::MODE)
      {
         details::CorrectSolutionFunctor<CouplingIndexType::MODE> func(res, cinfo, corr, matIdx);
         return func.corr;
      }
      // There is a single matrix
      else if(cinfo.indexType() == CouplingIndexType::SINGLE)
      {
         details::CorrectSolutionFunctor<CouplingIndexType::SINGLE> func(res, cinfo, corr, matIdx);
         return func.corr;
      }
      else
      {
         throw std::logic_error("Unknown coupling");
      }
   }

} // Equations
} // QuICC
