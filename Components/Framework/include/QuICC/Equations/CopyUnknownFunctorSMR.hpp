/**
 * @file CopyUnknownFunctorSMR.hpp
 * @brief Implementation of the CopyUnknown functor fpr SLOWEST_MULTI_RHS
 */

#ifndef QUICC_EQUATIONS_COPYUNKNOWNFUNCTORSMR_HPP
#define QUICC_EQUATIONS_COPYUNKNOWNFUNCTORSMR_HPP

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
#include "Arithmetics/Basic.hpp"
#include "QuICC/Equations/CopyUnknownFunctor.hpp"

namespace QuICC {

namespace Equations {

template <>
template <bool IsSet, typename TData, typename TField> void CopyUnknownFunctor<CouplingIndexType::SLOWEST_MULTI_RHS>::apply(const TField& field, TData& storage, const int start)
{
   const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);

   int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx);

   //Safety assertion
   assert(start >= 0);

   // Copy data
   for(int j = zeroCol; j < cols; j++)
   {
      const int rows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);
      for(int i = zeroRow; i < rows; i++)
      {
         if constexpr(IsSet)
         {
            // Copy field value into storage
            Arithmetics::setScalar(storage, i - zeroRow + start, j - zeroCol, field.comp(compId).point(i,j,matIdx));
         }
         else
         {
            // Add field value to storage
            Arithmetics::addScalar(storage, i - zeroRow + start, j - zeroCol, field.comp(compId).point(i,j,matIdx));
         }
      }
   }
}

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_COPYUNKNOWNFUNCTORSMR_HPP
