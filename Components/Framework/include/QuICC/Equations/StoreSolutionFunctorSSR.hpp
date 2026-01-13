/**
 * @file StoreSolutionFunctorSSR.hpp
 * @brief Implementation of the StoreSolution functor for SLOWEST_SINGLE_RHS
 */

#ifndef QUICC_EQUATIONS_STORESOLUTIONFUNCTORSSR_HPP
#define QUICC_EQUATIONS_STORESOLUTIONFUNCTORSSR_HPP

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
#include "Arithmetics/Basic.hpp"
#include "QuICC/Equations/IFieldEquation_decl.hpp"
#include "QuICC/Equations/StoreSolutionFunctor.hpp"

namespace QuICC {

namespace Equations {

template <>
template <typename TData, typename TField> void StoreSolutionFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS>::apply(TField& field, const TData& storage, const int start)
{
   int solStart;
   auto solution = init(solStart, storage, start, eq->couplingInfo(compId));

   const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
   int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx);

#if defined QUICC_MPI && defined QUICC_MPISPSOLVE
   // Add source data
   const auto& sRes = eq->res().sim();
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
   for(int j = 0; j < cols; j++)
   {
      j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
      if(corrDim > 0)
      {
         j_ -= corrDim;
      }
      const int rows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);
      for(int i = 0; i < rows; i++)
      {
         // Compute correct position
         l = start + j_ + i;

         // Copy timestep output into field
         MHDVariant dataPoint = Arithmetics::getScalar(*solution.ptr, l);
         dataPoint = eq->updateStoredSolution(dataPoint, compId, i, j, matIdx);
         field.rComp(compId).setPoint(dataPoint,i,j,matIdx);
      }
   }
#else
   // Copy data
   int k = solStart;
   for(int j = 0; j < cols; j++)
   {
      // Effective rows in case of non-uniform truncation
      int usedRows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);

      for(int i = 0; i < usedRows; i++)
      {
         // Copy timestep output into field
         MHDVariant dataPoint = Arithmetics::getScalar(*solution.ptr, k);
         dataPoint = eq->updateStoredSolution(dataPoint, compId, i, j, matIdx);
         field.rComp(compId).setPoint(dataPoint,i,j,matIdx);

         // increase linear storage counter
         k++;
      }
   }
#endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE
}

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_STORESOLUTIONFUNCTORSSR_HPP
