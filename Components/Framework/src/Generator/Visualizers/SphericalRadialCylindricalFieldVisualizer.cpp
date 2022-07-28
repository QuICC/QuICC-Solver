/**
 * @file SphericalRadialCylindricalFieldVisualizer.cpp
 * @brief Source of the implementation of the spherical vertical component field visualizer
 * @author Nicolò Lardelli \<nicolo.lardelli@erdw.ethz.ch\>
 */

// Configuration includes
//

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Generator/Visualizers/SphericalRadialCylindricalFieldVisualizer.hpp"

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/Generator/Visualizers/Kernels/Spherical/CylindricalRadialField.hpp"

namespace QuICC {

namespace Equations {

   SphericalRadialCylindricalFieldVisualizer::SphericalRadialCylindricalFieldVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IScalarEquation(spEqParams,spScheme,spBackend), mFieldType(FieldType::VECTOR)
   {
   }

   SphericalRadialCylindricalFieldVisualizer::~SphericalRadialCylindricalFieldVisualizer()
   {
   }

   void SphericalRadialCylindricalFieldVisualizer::setIdentity(const std::size_t vertName, const std::size_t fieldName)
   {
      // Set the name
      this->setName(vertName);

      // Store field name
      this->mFieldName = fieldName;

      // Set the variable requirements
      this->setRequirements();
   }

   void SphericalRadialCylindricalFieldVisualizer::setFieldType(const FieldType::Id type)
   {
      if(type == FieldType::GRADIENT)
      {
         throw std::logic_error("Z Component of gradient not implemented yet!");
      }

      this->mFieldType = type;
   }

   void SphericalRadialCylindricalFieldVisualizer::setCoupling()
   {
      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Nonlinear) = true;

      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::WRAPPER, 0, features);
   }

   void SphericalRadialCylindricalFieldVisualizer::initNLKernel(const bool force)
   {
      // Initialize if empty or forced
      if(force || !this->mspNLKernel)
      {
         // Initialize the physical kernel
         auto spNLKernel = std::make_shared<Physical::Kernel::Spherical::CylindricalRadialField>();
         spNLKernel->setVector(this->mFieldName, this->spVector(this->mFieldName));
         spNLKernel->setScalar(this->name(), this->spScalar(this->name()));
         spNLKernel->init(this->mFieldType, 1.0);
         this->mspNLKernel = spNLKernel;
      }
   }

   void SphericalRadialCylindricalFieldVisualizer::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::After::id());

      // Forward transform generates field values
      this->setForwardPathsType(FWD_IS_FIELD);

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add unknown to requirements: is scalar?
      auto& req = this->mRequirements.addField(this->name(), FieldRequirement(true, ss.spectral(), ss.physical()));
      req.enableSpectral();
      req.enablePhysical();

      // Add base vector field to requirements: is scalar?, need spectral?, need physical?, need diff?(, needCurl)
      auto& reqField = this->mRequirements.addField(this->mFieldName, FieldRequirement(false, ss.spectral(), ss.physical()));
      reqField.enableSpectral();
      if(this->mFieldType == FieldType::VECTOR) reqField.enablePhysical();
      if(this->mFieldType == FieldType::GRADIENT) reqField.enableGradient();
      if(this->mFieldType == FieldType::CURL) reqField.enableCurl();
   }

}
}
