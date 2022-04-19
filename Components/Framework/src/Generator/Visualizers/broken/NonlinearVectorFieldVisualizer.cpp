/**
 * @file NonlinearVectorFieldVisualizer.cpp
 * @brief Source of the implementation of a nonlinear vector field visualizer
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/Visualizers/NonlinearVectorFieldVisualizer.hpp"

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/PhysicalOperators/Cross.hpp"

namespace QuICC {

namespace Equations {

   NonlinearVectorFieldVisualizer::NonlinearVectorFieldVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme)
      : IVectorEquation(spEqParams,spScheme)
   {
   }

   NonlinearVectorFieldVisualizer::~NonlinearVectorFieldVisualizer()
   {
   }

   void NonlinearVectorFieldVisualizer::setIdentity(const std::size_t name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void NonlinearVectorFieldVisualizer::setNonlinearType(const NonlinearVectorVisualizerIds::Id type)
   {
      this->mNonlinearType = type;
   }

   void NonlinearVectorFieldVisualizer::setCoupling()
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

   void NonlinearVectorFieldVisualizer::computeNonlinear(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id compId) const
   {
      if(this->mNonlinearType == NonlinearVectorVisualizerIds::CYLINDER_TORPOL_ADVECTION)
      {
         switch(compId)
         {
            case(FieldComponents::Physical::R):
               Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::Z>::set(rNLComp, this->vector(this->name()).dom(0).curl(), this->vector(this->name()).dom(0).phys(), 1.0);
               break;
            case(FieldComponents::Physical::THETA):
               Physical::Cross<FieldComponents::Physical::Z,FieldComponents::Physical::R>::set(rNLComp, this->vector(this->name()).dom(0).curl(), this->vector(this->name()).dom(0).phys(), 1.0);
               break;
            case(FieldComponents::Physical::Z):
               Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::set(rNLComp, this->vector(this->name()).dom(0).curl(), this->vector(this->name()).dom(0).phys(), 1.0);
               break;
            default:
               assert(false);
               break;
         }
      }
   }

   void NonlinearVectorFieldVisualizer::useNonlinear(const Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id compId)
   {
      this->rUnknown().rDom(0).rPhys().rComp(compId).rData() = rNLComp.data();
   }

   void NonlinearVectorFieldVisualizer::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::After::id());

      if(this->mNonlinearType == NonlinearVectorVisualizerIds::CYLINDER_TORPOL_ADVECTION)
      {
         // Get reference to spatial scheme
         const auto& ss = this->ss();

         // Add unknown to requirements: is scalar?
         auto& req = this->mRequirements.addField(this->name(), FieldRequirement(false, ss.spectral(), ss.physical()));
         req.enableSpectral();
         req.enablePhysical();
         req.enableCurl();
      }
   }

   void NonlinearVectorFieldVisualizer::setNLComponents()
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
