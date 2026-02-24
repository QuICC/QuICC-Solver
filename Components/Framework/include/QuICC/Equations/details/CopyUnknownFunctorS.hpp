/**
 * @file CopyUnknownFunctorS.hpp
 * @brief Implementation of the CopyUnknown functor for SINGLE
 */

#ifndef QUICC_EQUATIONS_DETAILS_COPYUNKNOWNFUNCTORS_HPP
#define QUICC_EQUATIONS_DETAILS_COPYUNKNOWNFUNCTORS_HPP

// System includes
//
#include <memory>
#include <stdexcept>
#include <vector>

// Project includes
//
#include "Arithmetics/Basic.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/Equations/details/CopyUnknownFunctor.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <>
template <bool IsSet, typename TData, typename TField>
void CopyUnknownFunctor<CouplingIndexType::SINGLE>::apply(const TField& field,
   TData& storage, const int start)
{
   const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
   const auto& sRes = eq->res().sim();

   assert(matIdx == 0);

   // Safety assertion
   assert(start >= 0);

   // Copy data
   int l, k_, j_, dimK, dimJ;

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
            l = start + k_ + j_ + i;

            if constexpr (IsSet)
            {
               // Copy field value into storage
               Arithmetics::assignScalar<Arithmetics::Operation::Set>(storage,
                  l, field.comp(compId).point(i, j, k));
            }
            else
            {
               // Add field value to storage
               Arithmetics::assignScalar<Arithmetics::Operation::Plus>(storage,
                  l, field.comp(compId).point(i, j, k));
            }
         }
      }
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_COPYUNKNOWNFUNCTORS_HPP
