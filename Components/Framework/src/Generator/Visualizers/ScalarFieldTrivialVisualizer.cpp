/**
 * @file ScalarFieldTrivialVisualizer.cpp
 * @brief Source of the implementation of the basic scalar field visualizer with trivial solver
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/Visualizers/ScalarFieldTrivialVisualizer.hpp"

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Math.hpp"
#include "QuICC/SolveTiming/After.hpp"

namespace QuICC {

namespace Equations {

   ScalarFieldTrivialVisualizer::ScalarFieldTrivialVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IScalarEquation(spEqParams,spScheme,spBackend), mViewField(true), mViewGradient(false), mViewGradient2(false)
   {
   }

   ScalarFieldTrivialVisualizer::~ScalarFieldTrivialVisualizer()
   {
   }

   void ScalarFieldTrivialVisualizer::setIdentity(const std::size_t name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void ScalarFieldTrivialVisualizer::setFields(const bool viewField, const bool viewGradient, const bool viewGradient2)
   {
      this->mViewField = viewField;

      this->mViewGradient = viewGradient;

      this->mViewGradient2 = viewGradient2;
   }

   void ScalarFieldTrivialVisualizer::setCoupling()
   {
      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Nonlinear) = true;

      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, features);
   }

   void ScalarFieldTrivialVisualizer::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::After::id());

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
