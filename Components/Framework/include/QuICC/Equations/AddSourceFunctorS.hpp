/**
 * @file AddSourceFunctorS.hpp
 * @brief AddSourceFunctor implementation For SINGLE
 */

#ifndef QUICC_EQUATIONS_ADDSOURCEFUNCTORS_HPP
#define QUICC_EQUATIONS_ADDSOURCEFUNCTORS_HPP

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
   template <typename TData, typename TField> void AddSourceFunctor<CouplingIndexType::SINGLE>::apply(const TField& field, TData& storage, const int start)
   {
      const auto& info = eq->couplingInfo(compId);
      // Add source term if required
      if(info.hasSource())
      {
         const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
         const auto& sRes = eq->res().sim();
         assert(matIdx == 0);

         //int zeroRow = info.galerkinShift(matIdx,0);
         //int zeroCol = info.galerkinShift(matIdx,1);
         //int zeroBlock = info.galerkinShift(matIdx,2);

         //Safety assertion
         assert(start >= 0);

         // Add source term
         int l, k_, j_, dimK, dimJ;

         switch(sRes.ss().dimension())
         {
            case 3:
               dimK = sRes.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)*sRes.dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
               dimJ = sRes.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
               break;
            case 2:
               dimK = 1;
               dimJ = sRes.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
               break;
            case 1:
               dimK = 1;
               dimJ = 1;
               break;
            default:
               dimK = -1;
               dimJ = -1;
               throw  std::logic_error("Spatial scheme has unknown dimension!");
         }

         for(int k = 0; k < tRes.template dim<Dimensions::Data::DAT3D>(); k++)
         {
            k_ = tRes.template idx<Dimensions::Data::DAT3D>(k)*dimK;
            for(int j = 0; j < tRes.template dim<Dimensions::Data::DAT2D>(k); j++)
            {
               j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,k)*dimJ;
               for(int i = 0; i < sRes.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL); i++)
               {
                  // Compute correct position
                  l = start + k_ + j_ + i;

                  // Add source term
                  Arithmetics::addScalar(storage, l, eq->sourceTerm(compId, i, j, k));
               }
            }
         }
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_ADDSOURCEFUNCTORS_HPP
