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
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    * @param isSet   Set field? (or add)
    */
   template <typename TField, typename TData> void copyNonlinear(const Resolution& res, const CouplingInformation& cinfo, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const bool isSet = false);
   template <typename TField, typename TOperator, typename TData> void copyNonlinear(const Resolution& res, const CouplingInformation& cinfo, const TField& field, FieldComponents::Spectral::Id compId, const TOperator& op, TData& storage, const int matIdx, const int start, const bool isSet = false);

   template <typename TField, typename TOperator, typename TData> void copyNonlinear(const Resolution& res, const CouplingInformation& cinfo, const TField& field, FieldComponents::Spectral::Id compId, const TOperator& op, TData& storage, const int matIdx, const int start, const bool isSet)
   {
      assert((!cinfo.isGalerkin() || cinfo.indexType() != CouplingIndexType::SINGLE) && "Current version does not support galerkin basis");

      if constexpr(Arithmetics::is_view<TData>::value)
      {
         throw std::logic_error("Not yet implemented for View data");
      }
      else
      {
         // Temporary storage is required
         TData tmp;
         tmp = TData(cinfo.tauN(matIdx), cinfo.rhsCols(matIdx));

         // simply copy values from unknown
         copyUnknown(res, cinfo, field, compId, tmp, matIdx, 0, false, true, true);

         // Multiply nonlinear term by quasi-inverse
         applyQuasiInverse(op, compId, storage, start, matIdx, 0, tmp, isSet);
      }
   }

   template <typename TField, typename TData> void copyNonlinear(const Resolution& res, const CouplingInformation& cinfo, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const bool isSet)
   {
      assert((!cinfo.isGalerkin() || cinfo.indexType() != CouplingIndexType::SINGLE) && "Current version does not support galerkin basis");

      /// Nonlinear computation took place but no quasi-inverse is required
      if(cinfo.hasNonlinear())
      {
         // simply copy values from unknown
         copyUnknown(res, cinfo, field, compId, storage, matIdx, start, true, isSet, false);
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_COPYNONLINEAR_HPP
