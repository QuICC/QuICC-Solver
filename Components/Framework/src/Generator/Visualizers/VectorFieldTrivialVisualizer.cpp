/**
 * @file VectorFieldTrivialVisualizer.cpp
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
#include "QuICC/Generator/Visualizers/VectorFieldTrivialVisualizer.hpp"

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/Transform/Path/TorPol.hpp"

namespace QuICC {

namespace Equations {

   VectorFieldTrivialVisualizer::VectorFieldTrivialVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IVectorEquation(spEqParams,spScheme,spBackend), mViewField(true), mViewGradient(false), mViewCurl(false)
   {
   }

   VectorFieldTrivialVisualizer::~VectorFieldTrivialVisualizer()
   {
   }

   void VectorFieldTrivialVisualizer::setIdentity(const std::size_t name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void VectorFieldTrivialVisualizer::setFields(const bool viewField, const bool viewGradient, const bool viewCurl)
   {
      this->mViewField = viewField;

      this->mViewGradient = viewGradient;

      this->mViewCurl = viewCurl;
   }

   void VectorFieldTrivialVisualizer::setCoupling()
   {
      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Nonlinear) = true;

      if(this->ss().spectral().ONE() != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(this->ss().spectral().ONE(), CouplingInformation::TRIVIAL, 0, features);
      }

      if(this->ss().spectral().TWO() != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(this->ss().spectral().TWO(), CouplingInformation::TRIVIAL, 0, features);
      }

      if(this->ss().spectral().THREE() != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(this->ss().spectral().THREE(), CouplingInformation::TRIVIAL, 0, features);
      }
   }

   void VectorFieldTrivialVisualizer::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::After::id());

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add temperature to requirements: is scalar?
      auto& req = this->mRequirements.addField(this->name(), FieldRequirement(false, ss.spectral(), ss.physical()));
      req.enableSpectral();
      if(this->mViewField) req.enablePhysical();
      if(this->mViewGradient) req.enableGradient();
      if(this->mViewCurl) req.enableCurl();
   }

   void VectorFieldTrivialVisualizer::setNLComponents()
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
