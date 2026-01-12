/**
 * @file SetZeroNonlinear.hpp
 * @brief SetZeroNonlinear implementation
 */

#ifndef QUICC_EQUATIONS_SETZERONONLINEAR_HPP
#define QUICC_EQUATIONS_SETZERONONLINEAR_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "Arithmetics/Basic.hpp"
#include "QuICC/Equations/SetZeroNonlinearFunctor.hpp"
#include "QuICC/Equations/SetZeroNonlinearFunctorSSR.hpp"
#include "QuICC/Equations/SetZeroNonlinearFunctorSMR.hpp"
#include "QuICC/Equations/SetZeroNonlinearFunctorM.hpp"
#include "QuICC/Equations/SetZeroNonlinearFunctorS.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Set nonlinear spectral values to zero
    *
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData, typename TField> void setZeroNonlinear(const IFieldEquation& eq, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);


   template <typename TData, typename TField> void setZeroNonlinear(const IFieldEquation& eq, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // matIdx is the index of the slowest varying direction with a single RHS
      if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         SetZeroNonlinearFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS> func(eq, compId, matIdx);
         func.apply(field, storage, start);
      }
      // matIdx is the index of the slowest varying direction with multiple RHS
      else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
      {
         SetZeroNonlinearFunctor<CouplingIndexType::SLOWEST_MULTI_RHS> func(eq, compId, matIdx);
         func.apply(field, storage, start);
      }
      // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
      else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::MODE)
      {
         SetZeroNonlinearFunctor<CouplingIndexType::MODE> func(eq, compId, matIdx);
         func.apply(field, storage, start);
      }
      else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SINGLE)
      {
         SetZeroNonlinearFunctor<CouplingIndexType::SINGLE> func(eq, compId, matIdx);
         func.apply(field, storage, start);
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_SETZERONONLINEAR_HPP
