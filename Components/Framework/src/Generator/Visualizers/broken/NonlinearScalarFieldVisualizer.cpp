/**
 * @file NonlinearScalarFieldVisualizer.cpp
 * @brief Source of the implementation of a nonlinearscalar field visualizer
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/Visualizers/NonlinearScalarFieldVisualizer.hpp"

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Constants.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/TypeSelectors/TransformSelector.hpp"
#include "QuICC/PhysicalOperators/VelocityHeatAdvection.hpp"

namespace QuICC {

namespace Equations {

   NonlinearScalarFieldVisualizer::NonlinearScalarFieldVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme)
      : IScalarEquation(spEqParams,spScheme)
   {
   }

   NonlinearScalarFieldVisualizer::~NonlinearScalarFieldVisualizer()
   {
   }

   void NonlinearScalarFieldVisualizer::setIdentity(const std::size_t name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void NonlinearScalarFieldVisualizer::setNonlinearType(const NonlinearScalarVisualizerIds::Id type)
   {
      this->mNonlinearType = type;
   }

   void NonlinearScalarFieldVisualizer::setCoupling()
   {
      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Nonlinear) = true;

      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::WRAPPER, 0, features);
   }

   void NonlinearScalarFieldVisualizer::computeNonlinear(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id compId) const
   {
      if(this->mNonlinearType == NonlinearScalarVisualizerIds::CYLINDER_HEAT_ADVECTION)
      {
         // Assert on scalar component is used
         assert(compId == FieldComponents::Physical::SCALAR);
         Physical::VelocityHeatAdvection<FieldComponents::Physical::R,FieldComponents::Physical::THETA,FieldComponents::Physical::Z>::set(rNLComp, this->vector(PhysicalNames::Velocity::id()).dom(0).phys(), this->scalar(this->name()).dom(0).grad(), 1.0);
      }
   }

   void NonlinearScalarFieldVisualizer::useNonlinear(const Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id compId)
   {
      this->rUnknown().rDom(0).rPhys().rData() = rNLComp.data();
   }

   void NonlinearScalarFieldVisualizer::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::After::id());

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?, need curl?, need diff2?
      if(this->mNonlinearType == NonlinearScalarVisualizerIds::CYLINDER_HEAT_ADVECTION)
      {
         // Get reference to spatial scheme
         const auto& ss = this->ss();

         // Add temperature to requirements: is scalar?
         auto& req = this->mRequirements.addField(this->name(), FieldRequirement(true, ss.spectral(), ss.physical()));
         req.enableSpectral();
         req.enableGradient();

         // Add X velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
         auto& reqVel = this->mRequirements.addField(PhysicalNames::Velocity::id(), FieldRequirement(false, ss.spectral(), ss.physical()));
         reqVel.enableSpectral();
         reqVel.enablePhysical();
      }
   }

}
}
