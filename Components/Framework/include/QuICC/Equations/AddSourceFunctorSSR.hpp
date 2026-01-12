/**
 * @file AddSourceFunctorSMR.hpp
 * @brief AddSourceFunctor implementation For SLOWEST_SINGLE_RHS
 */

#ifndef QUICC_EQUATIONS_ADDSOURCEFUNCTORSSR_HPP
#define QUICC_EQUATIONS_ADDSOURCEFUNCTORSSR_HPP

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
   template <typename TData, typename TField> void AddSourceFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS>::apply(const TField& field, TData& storage, const int start)
   {
      const auto& info = eq->couplingInfo(compId);
      // Add source term if required
      if(info.hasSource())
      {
         const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
         const auto& sRes = eq->res().sim();
         int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx);
         int zeroRow = info.galerkinShift(matIdx,0);
         int zeroCol;
         if(sRes.ss().has(SpatialScheme::Feature::SpectralOrdering132))
         {
            zeroCol = info.galerkinShift(matIdx,2);
         } else
         {
            zeroCol = info.galerkinShift(matIdx,1);
         }

         //Safety assertion
         assert(start >= 0);

#if defined QUICC_MPI && defined QUICC_MPISPSOLVE
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
         for(int j = zeroCol; j < cols; j++)
         {
            j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
            if(corrDim > 0)
            {
               j_ -= corrDim;
            }
            const int rows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);
            for(int i = zeroRow; i < rows; i++)
            {
               // Compute correct position
               l = start + j_ + i;

               // Add source term
               Arithmetics::addScalar(storage, l, eq->sourceTerm(compId, i, j, matIdx));
            }
         }
#else
         // Add source term
         int k = start;
         for(int j = zeroCol; j < cols; j++)
         {
            // Effective rows in case of non-uniform truncation
            int usedRows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);

            for(int i = zeroRow; i < usedRows; i++)
            {
               // Add source term
               Arithmetics::addScalar(storage, k, eq->sourceTerm(compId, i, j, matIdx));

               // increase storage counter
               k++;
            }
         }
#endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_ADDSOURCEFUNCTORSSR_HPP
