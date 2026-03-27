/**
 * @file StoreSolutionFunctorM.hpp
 * @brief Implementation of the StoreSolution functor for MODE
 */

#ifndef QUICC_EQUATIONS_DETAILS_STORESOLUTIONFUNCTORM_HPP
#define QUICC_EQUATIONS_DETAILS_STORESOLUTIONFUNCTORM_HPP

// System includes
//


// Project includes
//
#include "Arithmetics/Basic.hpp"
#include "QuICC/Equations/details/StoreSolutionFunctor.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <>
template <typename TData, typename TField>
void StoreSolutionFunctor<CouplingIndexType::MODE>::apply(TField& field,
   const TData& storage, const int start)
{
   int solStart;
   auto solution = init(solStart, storage, start);

   const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
   // Get mode indexes
   ArrayI mode = tRes.mode(matIdx);
   int rows = field.comp(compId).slice(mode(0)).rows();

   // Copy data
   int k = solStart;
   for (int i = 0; i < rows; i++)
   {
      // Copy timestep output into field
      MHDVariant dataPoint = Arithmetics::getScalar(*solution.ptr, k);
      dataPoint =
         (*spUp)(dataPoint, i, mode(1), mode(0));
      field.rComp(compId).setPoint(dataPoint, i, mode(1), mode(0));

      // increase linear storage counter
      k++;
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_STORESOLUTIONFUNCTORM_HPP
