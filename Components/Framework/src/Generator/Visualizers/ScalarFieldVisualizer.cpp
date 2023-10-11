/**
 * @file ScalarFieldVisualizer.cpp
 * @brief Source of the implementation of the basic scalar field visualizer
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/Visualizers/ScalarFieldVisualizer.hpp"

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Math.hpp"
#include "QuICC/SolveTiming/After.hpp"

namespace QuICC {

namespace Equations {

   ScalarFieldVisualizer::ScalarFieldVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IScalarEquation(spEqParams,spScheme,spBackend), mViewField(true), mViewGradient(false), mViewGradient2(false)
   {
   }

   ScalarFieldVisualizer::~ScalarFieldVisualizer()
   {
   }

   void ScalarFieldVisualizer::setIdentity(const std::size_t name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void ScalarFieldVisualizer::setFields(const bool viewField, const bool viewGradient, const bool viewGradient2)
   {
      this->mViewField = viewField;

      this->mViewGradient = viewGradient;

      this->mViewGradient2 = viewGradient2;
   }

   void ScalarFieldVisualizer::setCoupling()
   {
      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Nonlinear) = true;

      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::WRAPPER, 0, features);
   }

   void ScalarFieldVisualizer::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::After::id());

      // Forward transform generates field values
      this->setForwardPathsType(FWD_IS_FIELD);

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add temperature to requirements: is scalar?
      auto& req = this->mRequirements.addField(this->name(), FieldRequirement(true, ss.spectral(), ss.physical()));
      req.enableSpectral();
      if(this->mViewField) req.enablePhysical();
      if(this->mViewGradient) req.enableGradient();
      if(this->mViewGradient2) req.enableGradient2();
   }

}
}
