/**
 * @file AddSource.hpp
 * @brief AddSource implementation
 */

#ifndef QUICC_EQUATIONS_ADDSOURCE_HPP
#define QUICC_EQUATIONS_ADDSOURCE_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/Equations/AddSourceFunctor.hpp"
#include "QuICC/Equations/AddSourceFunctorSMR.hpp"
#include "QuICC/Equations/AddSourceFunctorSSR.hpp"
#include "QuICC/Equations/AddSourceFunctorM.hpp"
#include "QuICC/Equations/AddSourceFunctorS.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Add source term
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData, typename TField> void addSource(const IFieldEquation& eq, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   template <typename TData, typename TField> void addSource(const IFieldEquation& eq, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // matIdx is the index of the slowest varying direction with a single RHS
      if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         AddSourceFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS> func(eq, compId, matIdx);
         func.apply(field, storage, start);
      }
      // matIdx is the index of the slowest varying direction with multiple RHS
      else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
      {
         AddSourceFunctor<CouplingIndexType::SLOWEST_MULTI_RHS> func(eq, compId, matIdx);
         func.apply(field, storage, start);
      }
      // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
      else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::MODE)
      {
         AddSourceFunctor<CouplingIndexType::MODE> func(eq, compId, matIdx);
         func.apply(field, storage, start);
      }
      // There is a single matrix
      else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SINGLE)
      {
         AddSourceFunctor<CouplingIndexType::SINGLE> func(eq, compId, matIdx);
         func.apply(field, storage, start);
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_ADDSOURCE_HPP
