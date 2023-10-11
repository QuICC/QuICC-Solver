/**
 * @file CartesianExactScalarState.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution in cartesian geometries
 */

// System includes
//

// Project includes
//
#include "QuICC/Generator/States/CartesianExactScalarState.hpp"
#include "Types/Typedefs.hpp"
#include "Types/Math.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/PhysicalKernels/DoNothing.hpp"
#include "QuICC/PhysicalKernels/MakeConstant.hpp"
#include "QuICC/PhysicalKernels/MakeRandom.hpp"
#include "QuICC/SpectralKernels/Set3DModes.hpp"

namespace QuICC {

namespace Equations {

   CartesianExactScalarState::CartesianExactScalarState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IScalarEquation(spEqParams,spScheme,spBackend)
   {
   }

   void CartesianExactScalarState::setIdentity(const std::size_t name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void CartesianExactScalarState::setPhysicalKernel(Physical::Kernel::SharedIPhysicalKernel spKernel)
   {
      this->mspPhysKernel = spKernel;
   }

   void CartesianExactScalarState::setPhysicalNoise(const MHDFloat level)
   {
      auto spKernel = std::make_shared<Physical::Kernel::MakeRandom>();
      spKernel->init(level);
      this->mspPhysKernel = spKernel;
   }

   void CartesianExactScalarState::setPhysicalConstant(const MHDFloat value)
   {
      auto spKernel = std::make_shared<Physical::Kernel::MakeConstant>();
      spKernel->init(value);
      this->mspPhysKernel = spKernel;
   }

   void CartesianExactScalarState::setSpectralModes(const Spectral::Kernel::Complex3DMapType& modes)
   {
      auto spKernel = std::make_shared<Spectral::Kernel::Set3DModes>(this->ss().has(SpatialScheme::Feature::ComplexSpectrum));
      spKernel->init(modes);
      this->setSrcKernel(FieldComponents::Spectral::SCALAR, spKernel);
   }

   void CartesianExactScalarState::setCoupling()
   {
      bool hasNL = false;
      bool hasSource = false;
      if(this->mspPhysKernel)
      {
         hasNL = true;
      }

      if(this->mSrcKernel.size() > 0)
      {
         hasSource = true;
      }

      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Nonlinear) = hasNL;
      features.at(CouplingFeature::Source) = hasSource;
      features.at(CouplingFeature::AllowExplicit) = false;

      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, features);
   }

   void CartesianExactScalarState::initNLKernel(const bool force)
   {
      // Initialize if empty or forced
      if(force || !this->mspNLKernel)
      {
         if(this->mspPhysKernel)
         {
            this->mspNLKernel = this->mspPhysKernel;
         } else if(this->mSrcKernel.size() > 0)
         {
            // Pure spectral state is used
            auto spNLKernel = std::make_shared<Physical::Kernel::DoNothing>();
            this->mspNLKernel = spNLKernel;
         } else
         {
            throw std::logic_error("Unknown exact state");
         }

         // Set resolution
         this->mspNLKernel->setResolution(this->spRes());
      }
   }

   void CartesianExactScalarState::setRequirements()
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
   }

}
}
