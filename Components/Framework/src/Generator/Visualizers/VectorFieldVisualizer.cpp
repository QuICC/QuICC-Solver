/**
 * @file VectorFieldVisualizer.cpp
 * @brief Source of the implementation of the basic vector field visualizer
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/Visualizers/VectorFieldVisualizer.hpp"

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Math.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/Transform/Path/TorPol.hpp"

namespace QuICC {

namespace Equations {

   VectorFieldVisualizer::VectorFieldVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IVectorEquation(spEqParams,spScheme,spBackend), mViewField(true), mViewGradient(false), mViewCurl(false)
   {
   }

   VectorFieldVisualizer::~VectorFieldVisualizer()
   {
   }

   void VectorFieldVisualizer::setIdentity(const std::size_t name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void VectorFieldVisualizer::setFields(const bool viewField, const bool viewGradient, const bool viewCurl)
   {
      this->mViewField = viewField;

      this->mViewGradient = viewGradient;

      this->mViewCurl = viewCurl;
   }

   void VectorFieldVisualizer::setCoupling()
   {
      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Nonlinear) = true;

      if(this->ss().spectral().ONE() != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(this->ss().spectral().ONE(), CouplingInformation::WRAPPER, 0, features);
      }

      if(this->ss().spectral().TWO() != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(this->ss().spectral().TWO(), CouplingInformation::WRAPPER, 0, features);
      }

      if(this->ss().spectral().THREE() != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(this->ss().spectral().THREE(), CouplingInformation::WRAPPER, 0, features);
      }
   }

   void VectorFieldVisualizer::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::After::id());

      // Forward transform generates field values
      this->setForwardPathsType(FWD_IS_FIELD);

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add temperature to requirements: is scalar?
      auto& req = this->mRequirements.addField(this->name(), FieldRequirement(false, ss.spectral(), ss.physical()));
      req.enableSpectral();
      if(this->mViewField) req.enablePhysical();
      if(this->mViewGradient) req.enableGradient();
      if(this->mViewCurl) req.enableCurl();
   }

   void VectorFieldVisualizer::setNLComponents()
   {
      std::size_t pathId;
      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         pathId = Transform::Path::TorPol::id();
      }
      else
      {
         pathId = 0;
      }

      if(this->ss().spectral().ONE() != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLComponent(this->ss().spectral().ONE(), pathId);
      }

      if(this->ss().spectral().TWO() != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLComponent(this->ss().spectral().TWO(), pathId);
      }

      if(this->ss().spectral().THREE() != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLComponent(this->ss().spectral().THREE(), pathId);
      }
   }

}
}
