/**
 * @file IFieldEquation.cpp
 * @brief Source of generic field equation interface
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Equations/IFieldEquation.hpp"

// Project includes
//

namespace QuICC {

namespace Equations {

   IFieldEquation::IFieldEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IEquation(spEqParams, spScheme, spBackend)
   {
   }

   IFieldEquation::~IFieldEquation()
   {
   }

   void IFieldEquation::applyConstraint(FieldComponents::Spectral::Id compId)
   {
      // Use source kernel
      if(this->mConstraintKernel.count(compId) > 0)
      {
         return this->mConstraintKernel.find(compId)->second->apply();
      }
   }

   MHDVariant IFieldEquation::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
   {
      // Use source kernel
      if(this->mSrcKernel.count(compId) > 0)
      {
         return this->mSrcKernel.find(compId)->second->compute(i, j, k);
      } else
      {
         // This implementation should never get called!
         throw std::logic_error("Activated source term without implementation!");

         return MHDVariant();
      }
   }

   MHDVariant IFieldEquation::boundaryValue(FieldComponents::Spectral::Id, const int, const int, const int) const
   {
      // This implementation should never get called!
      throw std::logic_error("Activated boundary value without implementation!");

      return MHDVariant();
   }
}
}
