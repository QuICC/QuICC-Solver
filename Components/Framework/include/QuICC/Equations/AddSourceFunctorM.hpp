/**
 * @file AddSourceFunctorM.hpp
 * @brief AddSourceFunctor implementation For MODE
 */

#ifndef QUICC_EQUATIONS_ADDSOURCEFUNCTORM_HPP
#define QUICC_EQUATIONS_ADDSOURCEFUNCTORM_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/AddSourceFunctor.hpp"
#include "Arithmetics/Basic.hpp"

namespace QuICC {

namespace Equations {

   template <>
   template <typename TData, typename TField> void AddSourceFunctor<CouplingIndexType::MODE>::apply(const TField& field, TData& storage, const int start)
   {
      const auto& info = eq->couplingInfo(compId);
      // Add source term if required
      if(info.hasSource())
      {
         const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
         //Safety assertion
         assert(start >= 0);

         // Get mode indexes
         ArrayI mode = tRes.mode(matIdx);
         int rows = field.comp(compId).slice(mode(0)).rows();
         int zeroRow = info.galerkinShift(matIdx,0);

         // Copy data
         int k = start;
         for(int i = zeroRow; i < rows; i++)
         {
            // Add source term
            Arithmetics::addScalar(storage, k, eq->sourceTerm(compId, i, mode(1), mode(0)));

            // increase storage counter
            k++;
         }
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_ADDSOURCEFUNCTORM_HPP
