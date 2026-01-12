/**
 * @file CopyUnknownFunctorM.hpp
 * @brief Implementation of the CopyUnknown functor for MODE
 */

#ifndef QUICC_EQUATIONS_COPYUNKNOWNFUNCTORM_HPP
#define QUICC_EQUATIONS_COPYUNKNOWNFUNCTORM_HPP

// System includes
//
#include <vector>
#include <memory>
#include <stdexcept>


// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/Equations/SetZeroNonlinear.hpp"
#include "Arithmetics/Basic.hpp"
#include "QuICC/Equations/CopyUnknownFunctor.hpp"

namespace QuICC {

namespace Equations {

template <>
template <bool IsSet, typename TData, typename TField> void CopyUnknownFunctor<CouplingIndexType::MODE>::apply(const TField& field, TData& storage, const int start)
{
   const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
   const auto& sRes = eq->res().sim();

   //Safety assertion
   assert(start >= 0);

   // Get mode indexes
   ArrayI mode = tRes.mode(matIdx);
   int rows = field.comp(compId).slice(mode(0)).rows();

   bool isUsed;
   if(sRes.ss().has(SpatialScheme::Feature::FourierIndex23))
   {
      // Filter out complex conjugate modes to be safe
      isUsed = !(mode(3) == 0 && mode(2) > sRes.dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2);
   }
   else
   {
      isUsed = true;
   }

   if(isUsed)
   {
      // Copy data
      int k = start;
      for(int i = zeroRow; i < rows; i++)
      {
         if constexpr(IsSet)
         {
            // Copy field value into storage
            Arithmetics::setScalar(storage, k, field.comp(compId).point(i,mode(1),mode(0)));
         }
         else
         {
            // Add field value to storage
            Arithmetics::addScalar(storage, k, field.comp(compId).point(i,mode(1),mode(0)));
         }

         // increase storage counter
         k++;
      }
   }
   else
   {
      setZeroNonlinear(*eq, field, compId, storage, matIdx, start);
   }
}

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_COPYUNKNOWNFUNCTORM_HPP
