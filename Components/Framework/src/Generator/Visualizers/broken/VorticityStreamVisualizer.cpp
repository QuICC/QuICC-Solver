/**
 * @file VorticityStreamVisualizer.cpp
 * @brief Source of the implementation of the streamfunction to vorticity visualizer
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "QuICC/Generator/Visualizers/VorticityStreamVisualizer.hpp"

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/PhysicalNames/Streamfunction.hpp"
#include "QuICC/PhysicalNames/Vorticity.hpp"

namespace QuICC {

namespace Equations {

   VorticityStreamVisualizer::VorticityStreamVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme)
      : IScalarEquation(spEqParams,spScheme), mViewField(true), mViewGradient(false)
   {
   }

   VorticityStreamVisualizer::~VorticityStreamVisualizer()
   {
   }

   void VorticityStreamVisualizer::setFields(const bool viewField, const bool viewGradient)
   {
      this->mViewField = viewField;

      this->mViewGradient = viewGradient;

      // Set the variable requirements
      this->setRequirements();
   }

   void VorticityStreamVisualizer::setCoupling()
   {
      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::AllowExplicit) = false;

      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::WRAPPER, 0, features);
   }

   void VorticityStreamVisualizer::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::Vorticity::id());

      // Set solver timing
      this->setSolveTiming(SolveTiming::After::id());

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add temperature to requirements: is scalar?
      auto& reqVor = this->mRequirements.addField(PhysicalNames::Vorticity::id(), FieldRequirement(true, ss.spectral(), ss.physical()));
      reqVor.enableSpectral();
      if(this->mViewField) reqVor.enablePhysical();
      if(this->mViewGradient) reqVor.enableGradient();

      // Add streamfunction requirements: is scalar?
      auto& reqStr = this->mRequirements.addField(PhysicalNames::Streamfunction::id(), FieldRequirement(true, ss.spectral(), ss.physical()));
      reqStr.enableSpectral();
   }

}
}
