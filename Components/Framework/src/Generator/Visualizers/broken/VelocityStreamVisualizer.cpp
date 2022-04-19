/**
 * @file VelocityStreamVisualizer.cpp
 * @brief Source of the implementation of the streamfunction + axial velocity to velocity field visualizer
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/Visualizers/VelocityStreamVisualizer.hpp"

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/SolveTiming/After.hpp"

namespace QuICC {

namespace Equations {

   VelocityStreamVisualizer::VelocityStreamVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme)
      : IVectorEquation(spEqParams,spScheme), mViewField(true), mViewGradient(false)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   VelocityStreamVisualizer::~VelocityStreamVisualizer()
   {
   }

   void VelocityStreamVisualizer::setIdentity(const std::size_t name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void VelocityStreamVisualizer::setFields(const bool viewField, const bool viewGradient)
   {
      this->mViewField = viewField;

      this->mViewGradient = viewGradient;
   }

   void VelocityStreamVisualizer::setCoupling()
   {
      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::AllowExplicit) = false;

      SpectralComponent_range specRange = this->spectralRange();
      SpectralComponent_iterator specIt;
      for(specIt = specRange.first; specIt != specRange.second; ++specIt)
      {
         this->defineCoupling(*specIt, CouplingInformation::WRAPPER, 0, features);
      }
   }

   void VelocityStreamVisualizer::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::After::id());

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add temperature to requirements: is scalar?
      auto& req = this->mRequirements.addField(this->name(), FieldRequirement(false, ss.spectral(), ss.physical()));
      if(this->mViewField) req.enablePhysical();
      if(this->mViewGradient) req.enableGradient();
   }

   void VelocityStreamVisualizer::setNLComponents()
   {
      if(this->ss().spectral().ONE() != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLComponent(this->ss().spectral().ONE(), 0);
      }

      if(this->ss().spectral().TWO() != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLComponent(this->ss().spectral().TWO(), 0);
      }

      if(this->ss().spectral().THREE() != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLComponent(this->ss().spectral().THREE(), 0);
      }
   }

}
}
