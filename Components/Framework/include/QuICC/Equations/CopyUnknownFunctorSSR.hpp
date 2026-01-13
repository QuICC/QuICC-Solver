/**
 * @file CopyUnknownFunctorSSR.hpp
 * @brief Implementation of the CopyUnknown functor for SLOWEST_SINGLE_RHS
 */

#ifndef QUICC_EQUATIONS_COPYUNKNOWNFUNCTORSSR_HPP
#define QUICC_EQUATIONS_COPYUNKNOWNFUNCTORSSR_HPP

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
template <bool IsSet, typename TData, typename TField> void CopyUnknownFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS>::apply(const TField& field, TData& storage, const int start)
{
   const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);

   int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx);

   //Safety assertion
   assert(start >= 0);

#if defined QUICC_MPI && defined QUICC_MPISPSOLVE
   const auto& sRes = eq->res().sim();
   // Add source data
   int l;
   int j_;
   int dimI = sRes.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
   int corrDim;
   if((sRes.ss().has(SpatialScheme::Feature::ShellGeometry) || sRes.ss().has(SpatialScheme::Feature::SphereGeometry)) &&
         sRes.ss().has(SpatialScheme::Feature::SpectralOrdering123) &&
         sRes.ss().has(SpatialScheme::Feature::SpectralMatrix2D))
   {
      corrDim = tRes.template idx<Dimensions::Data::DAT3D>(matIdx)*dimI;
   }
   if constexpr(IsSet)
   {
      ///\mhdBug This is overkill
      // Set storage to zero
      setZeroNonlinear(eq, field, compId, storage, matIdx, start);
   }

   for(int j = zeroCol; j < cols; j++)
   {
      j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
      const int rows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);
      if(corrDim > 0)
      {
         j_ -= corrDim;
      }
      for(int i = zeroRow; i < rows; i++)
      {
         // Compute correct position
         l = start + j_ + i;

         if constexpr(IsSet)
         {
            // Copy field value into storage
            Arithmetics::setScalar(storage, l, field.comp(compId).point(i,j,matIdx));
         }
         else
         {
            // Copy field value into storage
            Arithmetics::addScalar(storage, l, field.comp(compId).point(i,j,matIdx));
         }
      }
   }
#else
   // Copy data
   int k = start;
   for(int j = zeroCol; j < cols; j++)
   {
      // Effective rows in case of non-uniform truncation
      int usedRows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);

      for(int i = zeroRow; i < usedRows; i++)
      {
         if constexpr(IsSet)
         {
            // Copy field value into storage
            Arithmetics::setScalar(storage, k, field.comp(compId).point(i,j,matIdx));
         }
         else
         {
            // Add field value to storage
            Arithmetics::addScalar(storage, k, field.comp(compId).point(i,j,matIdx));
         }

         // increase storage counter
         k++;
      }
   }
#endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE
}

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_COPYUNKNOWNFUNCTORSSR_HPP
