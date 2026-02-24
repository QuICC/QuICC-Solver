/**
 * @file CopyNonlinear.hpp
 * @brief Base for the implementation of a vector equation
 */

#ifndef QUICC_EQUATIONS_COPYNONLINEAR_HPP
#define QUICC_EQUATIONS_COPYNONLINEAR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CopyUnknown.hpp"
#include "QuICC/Equations/ApplyQuasiInverse.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Transfer nonlinear spectral values from unknown to solver
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    * @param isSet   Set field? (or add)
    */
   template <typename TEquation, typename TData> void copyNonlinear(const TEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const bool isSet = false);

   template <typename TEquation, typename TData> void copyNonlinear(const TEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const bool isSet)
   {
      const auto& info = eq.couplingInfo(compId);
      assert((!info.isGalerkin() || info.indexType() != CouplingIndexType::SINGLE) && "Current version does not support galerkin basis");

      // Check if a nonlinear computation took place and a quasi-inverse has to be applied
      if(info.hasNonlinear() && info.hasQuasiInverse())
      {
         if constexpr(Arithmetics::is_view<TData>::value)
         {
            throw std::logic_error("Not yet implemented for View data");
         }
         else
         {
            // Temporary storage is required
            TData tmp;
            tmp = TData(info.tauN(matIdx), info.rhsCols(matIdx));

            // simply copy values from unknown
            std::visit(
                  [&](auto&& p)
                  {
                  copyUnknown(eq, p->dom(0).perturbation(), compId, tmp, matIdx, 0, false, true, true);
                  }, eq.spUnknown());

            // Multiply nonlinear term by quasi-inverse
            applyQuasiInverse(eq, compId, storage, start, matIdx, 0, tmp, isSet);
         }
      }
      /// Nonlinear computation took place but no quasi-inverse is required
      else if(info.hasNonlinear())
      {
         // simply copy values from unknown
         std::visit(
               [&](auto&& p)
               {
                  copyUnknown(eq, p->dom(0).perturbation(), compId, storage, matIdx, start, true, isSet, false);
               }, eq.spUnknown());
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_COPYNONLINEAR_HPP
