/**
 * @file StoreSolutionFunctorS.hpp
 * @brief Implementation of the StoreSolution functor for SINGLE
 */

#ifndef QUICC_EQUATIONS_DETAILS_STORESOLUTIONFUNCTORS_HPP
#define QUICC_EQUATIONS_DETAILS_STORESOLUTIONFUNCTORS_HPP

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
void StoreSolutionFunctor<CouplingIndexType::SINGLE>::apply(TField& field,
   const TData& storage, const int start)
{
   int solStart;
   auto solution = init(solStart, storage, start);

   const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
   assert(matIdx == 0);

   // Copy data
   int l, k_, j_, dimK, dimJ;

   const auto& sRes = res.sim();
   switch (sRes.ss().dimension())
   {
   case 3:
      dimK =
         sRes.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL) *
         sRes.dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
      dimJ =
         sRes.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      break;
   case 2:
      dimK = 1;
      dimJ =
         sRes.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      break;
   case 1:
      dimK = 1;
      dimJ = 1;
      break;
   default:
      dimK = -1;
      dimJ = -1;
      throw std::logic_error("Spatial scheme has unknown dimension!");
   }

   for (int k = 0; k < tRes.template dim<Dimensions::Data::DAT3D>(); k++)
   {
      k_ = tRes.template idx<Dimensions::Data::DAT3D>(k) * dimK;
      for (int j = 0; j < tRes.template dim<Dimensions::Data::DAT2D>(k); j++)
      {
         j_ = tRes.template idx<Dimensions::Data::DAT2D>(j, k) * dimJ;
         for (int i = 0; i < sRes.dim(Dimensions::Simulation::SIM1D,
                                Dimensions::Space::SPECTRAL);
            i++)
         {
            // Compute correct position
            l = solStart + k_ + j_ + i;

            // Copy timestep output into field
            MHDVariant dataPoint = Arithmetics::getScalar(*solution.ptr, l);
            dataPoint = (*spUp)(dataPoint, i, j, k);
            field.rComp(compId).setPoint(dataPoint, i, j, k);
         }
      }
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_STORESOLUTIONFUNCTORS_HPP
