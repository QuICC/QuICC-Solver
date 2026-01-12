/**
 * @file SetZeroNonlinearFunctorS.hpp
 * @brief SetZeroNonlinear implementation
 */

#ifndef QUICC_EQUATIONS_SETZERONONLINEARFUNCTORS_HPP
#define QUICC_EQUATIONS_SETZERONONLINEARFUNCTORS_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/Equations/SetZeroNonlinearFunctor.hpp"

namespace QuICC {

namespace Equations {

template <>
   template <typename TData, typename TField> void SetZeroNonlinearFunctor<CouplingIndexType::SINGLE>::apply(const TField& field, TData& storage, const int start)
   {
      const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
      const auto& sRes = eq->res().sim();
      //Safety assertion
      assert(matIdx == 0);
      assert(start >= 0);

#if defined QUICC_MPI && defined QUICC_MPISPSOLVE
      for(int k = 0; k < eq->couplingInfo(compId).galerkinN(matIdx); ++k)
      {
         // Set field to zero
         Arithmetics::setZero(storage, k + start);
      }
#else
      // Set data to zero
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

               // Set field to zero
               Arithmetics::setZero(storage, l);
            }
         }
      }
#endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_SETZERONONLINEARFUNCTORS_HPP
