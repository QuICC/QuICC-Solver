/**
 * @file IFieldEquation.cpp
 * @brief Source of generic field equation interface
 */

// System includes
//

// Project includes
//
#include "QuICC/Equations/IFieldEquation.hpp"

namespace QuICC {

namespace Equations {

   IFieldEquation::IFieldEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IEquation(spEqParams, spScheme, spBackend)
   {
   }

   IFieldEquation::IFieldEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend, std::shared_ptr<EquationOptions> spOptions)
      : IEquation(spEqParams, spScheme, spBackend, spOptions)
   {
   }

   bool IFieldEquation::applyConstraint(FieldComponents::Spectral::Id compId)
   {
      bool changedSolution = false;

      // Use source kernel
      if(this->mConstraintKernel.count(compId) > 0)
      {
         changedSolution = true;
         this->mConstraintKernel.find(compId)->second->apply();
      }

      return changedSolution;
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
} // Equations
} // QuICC
