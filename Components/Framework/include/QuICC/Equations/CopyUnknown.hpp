/**
 * @file CopyUnknown.hpp
 * @brief CopyUnknown function
 */

#ifndef QUICC_EQUATIONS_COPYUNKNOWN_HPP
#define QUICC_EQUATIONS_COPYUNKNOWN_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/IEquation.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "Arithmetics/Basic.hpp"
#include "QuICC/Equations/CopyUnknownFunctor.hpp"
#include "QuICC/Equations/CopyUnknownFunctorSMR.hpp"
#include "QuICC/Equations/CopyUnknownFunctorSSR.hpp"
#include "QuICC/Equations/CopyUnknownFunctorS.hpp"
#include "QuICC/Equations/CopyUnknownFunctorM.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Copy unknown spectral values to solver
    *
    * @param eq         Equation to work on
    * @param compId     Component ID
    * @param storage    Storage for the equation values
    * @param matIdx     Index of the given data
    * @param start      Start index for the storage
    * @param useShift   Use galerkin shifts
    * @param isSet      Arithmetic operation is set
    */
   template <typename TData, typename TField> void copyUnknown(const IFieldEquation& eq, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const bool useShift, const bool isSet);

   template <typename TData, typename TField> void copyUnknown(const IFieldEquation& eq, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const bool useShift, const bool isSet)
   {
      // matIdx is the index of the slowest varying direction with a single RHS
      if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         CopyUnknownFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS> func(eq, compId, matIdx, useShift);

         // Copy data
         if(isSet)
         {
            func.apply<true>(field, storage, start);
         }
         else
         {
            func.apply<false>(field, storage, start);
         }
      }
      // matIdx is the index of the slowest varying direction with multiple RHS
      else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
      {
         CopyUnknownFunctor<CouplingIndexType::SLOWEST_MULTI_RHS> func(eq, compId, matIdx, useShift);

         // Copy data
         if(isSet)
         {
            func.apply<true>(field, storage, start);
         }
         else
         {
            func.apply<false>(field, storage, start);
         }
      }
      // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
      else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::MODE)
      {
         CopyUnknownFunctor<CouplingIndexType::MODE> func(eq, compId, matIdx, useShift);

         // Copy data
         if(isSet)
         {
            func.apply<true>(field, storage, start);
         }
         else
         {
            func.apply<false>(field, storage, start);
         }
      }
      // There is a single matrix
      else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SINGLE)
      {
         CopyUnknownFunctor<CouplingIndexType::SINGLE> func(eq, compId, matIdx, useShift);

         // Copy data
         if(isSet)
         {
            func.apply<true>(field, storage, start);
         }
         else
         {
            func.apply<false>(field, storage, start);
         }
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_COPYUNKNOWN_HPP
