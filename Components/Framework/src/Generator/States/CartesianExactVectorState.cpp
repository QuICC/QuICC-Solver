/**
 * @file CartesianExactVectorState.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution in cartesian geometries
 */

// System includes
//

// Project includes
//
#include "QuICC/Generator/States/CartesianExactVectorState.hpp"
#include "QuICC/Typedefs.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/PhysicalKernels/DoNothing.hpp"
#include "QuICC/PhysicalKernels/MakeConstant.hpp"
#include "QuICC/PhysicalKernels/MakeRandom.hpp"
#include "QuICC/SpectralKernels/Set3DModes.hpp"
#include "QuICC/Transform/Path/TorPol.hpp"

namespace QuICC {

namespace Equations {

   CartesianExactVectorState::CartesianExactVectorState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IVectorEquation(spEqParams,spScheme,spBackend)
   {
   }

   void CartesianExactVectorState::setIdentity(const std::size_t name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void CartesianExactVectorState::setPhysicalKernel(Physical::Kernel::SharedIPhysicalKernel spKernel)
   {
      this->mspPhysKernel = spKernel;
   }

   void CartesianExactVectorState::setPhysicalNoise(const MHDFloat level)
   {
      auto spKernel = std::make_shared<Physical::Kernel::MakeRandom>();
      spKernel->init(level);
      this->mspPhysKernel = spKernel;
   }

   void CartesianExactVectorState::setPhysicalConstant(const MHDFloat value)
   {
      auto spKernel = std::make_shared<Physical::Kernel::MakeConstant>();
      spKernel->init(value);
      this->mspPhysKernel = spKernel;
   }

   void CartesianExactVectorState::setSpectralModes(const FieldComponents::Spectral::Id compId, const Spectral::Kernel::Complex3DMapType& modes)
   {
      auto spKernel = std::make_shared<Spectral::Kernel::Set3DModes>(this->ss().has(SpatialScheme::Feature::ComplexSpectrum));
      spKernel->init(modes);
      this->setSrcKernel(compId, spKernel);
   }

   void CartesianExactVectorState::setCoupling()
   {
      bool hasNL = false;
      bool hasSource = false;
      if(this->mspPhysKernel)
      {
         hasNL = true;
      }

      if(this->mSrcKernel.size())
      {
         hasSource = true;
      }

      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Nonlinear) = true;
      features.at(CouplingFeature::AllowExplicit) = false;
      features.at(CouplingFeature::AllowExplicit) = false;

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

   void CartesianExactVectorState::initNLKernel(const bool force)
   {
      // Initialize if empty or forced
      if(force || !this->mspNLKernel)
      {
         if(this->mspPhysKernel)
         {
            this->mspNLKernel = this->mspPhysKernel;

         } else if(this->mSrcKernel.size() > 0)
         {
            // Pur spectral state is used
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

   void CartesianExactVectorState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::After::id());

      // Forward transform generates field values
      this->setForwardPathsType(FWD_IS_FIELD);

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add unknown to requirements: is scalar?
      auto& req = this->mRequirements.addField(this->name(), FieldRequirement(false, ss.spectral(), ss.physical()));
      req.enableSpectral();
      req.enablePhysical();
   }

   void CartesianExactVectorState::setNLComponents()
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
