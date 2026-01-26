/**
 * @file SetBoundaryValueFunctorSSR.hpp
 * @brief SetBoundaryValue implementation
 */

#ifndef QUICC_EQUATIONS_DETAILS_SETBOUNDARYVALUEFUNCTORSSR_HPP
#define QUICC_EQUATIONS_DETAILS_SETBOUNDARYVALUEFUNCTORSSR_HPP

// System includes
//

// Project includes
//
#include "Arithmetics/Basic.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/Equations/details/SetBoundaryValueFunctor.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <>
template <typename TData, typename TField>
void SetBoundaryValueFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS>::apply(
   const TField& field, TData& storage, const int start)
{
   const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
   const auto& sRes = eq->res().sim();
   const auto& info = eq->couplingInfo(compId);
   int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx);
   int zeroRow = info.galerkinShift(matIdx, 0);
   int zeroCol;
   if (sRes.ss().has(SpatialScheme::Feature::SpectralOrdering132))
   {
      zeroCol = info.galerkinShift(matIdx, 2);
   }
   else
   {
      zeroCol = info.galerkinShift(matIdx, 1);
   }

   // Safety assertion
   assert(start >= 0);

#if defined QUICC_MPI && defined QUICC_MPISPSOLVE
   // Set boundary value
   int l;
   int j_;
   int dimI =
      sRes.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
   int corrDim;
   if ((sRes.ss().has(SpatialScheme::Feature::ShellGeometry) ||
          sRes.ss().has(SpatialScheme::Feature::SphereGeometry)) &&
       sRes.ss().has(SpatialScheme::Feature::SpectralOrdering123) &&
       sRes.ss().has(SpatialScheme::Feature::SpectralMatrix2D))
   {
      corrDim = tRes.template idx<Dimensions::Data::DAT3D>(matIdx) * dimI;
   }
   for (int j = zeroCol; j < cols; j++)
   {
      j_ = tRes.template idx<Dimensions::Data::DAT2D>(j, matIdx) * dimI;
      if (corrDim > 0)
      {
         j_ -= corrDim;
      }
      const int rows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);
      for (int i = zeroRow; i < rows; i++)
      {
         // Compute correct position
         l = start + j_ + i;

         // Add source term
         Arithmetics::assignScalar<Arithmetics::Operation::Set>(storage, l,
            eq->boundaryValue(compId, i, j, matIdx));
      }
   }
#else
   // Set boundary value
   int k = start;
   for (int j = zeroCol; j < cols; j++)
   {
      // Effective rows in case of non-uniform truncation
      int usedRows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);

      for (int i = zeroRow; i < usedRows; i++)
      {
         // Add source term
         Arithmetics::assignScalar<Arithmetics::Operation::Set>(storage, k,
            eq->boundaryValue(compId, i, j, matIdx));

         // increase storage counter
         k++;
      }
   }
#endif // defined QUICC_MPI && defined QUICC_MPISPSOLVE
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_SETBOUNDARYVALUEFUNCTORSSR_HPP
