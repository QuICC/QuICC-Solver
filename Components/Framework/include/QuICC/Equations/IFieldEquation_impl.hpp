/**
 * @file IFieldEquation.hpp
 * @brief Base building block for the implementation of an equation
 */

#ifndef QUICC_EQUATIONS_IFIELDEQUATION_IMPL_HPP
#define QUICC_EQUATIONS_IFIELDEQUATION_IMPL_HPP

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
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "Arithmetics/Basic.hpp"
#include "QuICC/Equations/IFieldEquation_decl.hpp"
#include "QuICC/Equations/StoreSolutionFunctor.hpp"
#include "QuICC/Equations/StoreSolutionFunctorSSR.hpp"
#include "QuICC/Equations/StoreSolutionFunctorSMR.hpp"
#include "QuICC/Equations/StoreSolutionFunctorM.hpp"
#include "QuICC/Equations/StoreSolutionFunctorS.hpp"

namespace QuICC {

namespace Equations {

   inline MHDVariant IFieldEquation::updateStoredSolution(const MHDVariant newData, FieldComponents::Spectral::Id, const int, const int, const int)
   {
      return newData;
   }

   template<typename TData, typename TField>
      void IFieldEquation::storeSolutionImpl(TField& field, FieldComponents::Spectral::Id compId, const TData& storage, const int matIdx, const int start)
   {
      if(this->couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         StoreSolutionFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS> func(*this, compId, matIdx);
         func.apply(field, storage, start);
      }
      else if(this->couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
      {
         StoreSolutionFunctor<CouplingIndexType::SLOWEST_MULTI_RHS> func(*this, compId, matIdx);
         func.apply(field, storage, start);
      }
      else if(this->couplingInfo(compId).indexType() == CouplingIndexType::MODE)
      {
         StoreSolutionFunctor<CouplingIndexType::MODE> func(*this, compId, matIdx);
         func.apply(field, storage, start);
      }
      else if(this->couplingInfo(compId).indexType() == CouplingIndexType::SINGLE)
      {
         StoreSolutionFunctor<CouplingIndexType::SINGLE> func(*this, compId, matIdx);
         func.apply(field, storage, start);
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_IFIELDEQUATION_IMPL_HPP
