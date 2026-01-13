/**
 * @file StoreSolutionFunctorM.hpp
 * @brief Implementation of the StoreSolution functor for MODE
 */

#ifndef QUICC_EQUATIONS_STORESOLUTIONFUNCTORM_HPP
#define QUICC_EQUATIONS_STORESOLUTIONFUNCTORM_HPP

// System includes
//


// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "Arithmetics/Basic.hpp"
#include "QuICC/Equations/IFieldEquation_decl.hpp"
#include "QuICC/Equations/StoreSolutionFunctor.hpp"

namespace QuICC {

namespace Equations {

template <>
template <typename TData, typename TField> void StoreSolutionFunctor<CouplingIndexType::MODE>::apply(TField& field, const TData& storage, const int start)
{
   int solStart;
   auto solution = init(solStart, storage, start, eq->couplingInfo(compId));

   const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
   // Get mode indexes
   ArrayI mode = tRes.mode(matIdx);
   int rows = field.comp(compId).slice(mode(0)).rows();

   // Copy data
   int k = solStart;
   for(int i = 0; i < rows; i++)
   {
      // Copy timestep output into field
      MHDVariant dataPoint = Arithmetics::getScalar(*solution.ptr, k);
      dataPoint = eq->updateStoredSolution(dataPoint, compId, i, mode(1), mode(0));
      field.rComp(compId).setPoint(dataPoint,i,mode(1),mode(0));

      // increase linear storage counter
      k++;
   }
}

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_STORESOLUTIONFUNCTORM_HPP
