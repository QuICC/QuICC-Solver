/**
 * @file RandomVectorState.cpp
 * @brief Source of the implementation of the general random scalar state equation
 */

// Configuration includes
//

// System includes
//
#include <time.h>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Generator/States/RandomVectorState.hpp"

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Math.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/Transform/Path/TorPol.hpp"

namespace QuICC {

namespace Equations {

   RandomVectorState::RandomVectorState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IVectorEquation(spEqParams,spScheme,spBackend)
   {
   }

   RandomVectorState::~RandomVectorState()
   {
   }

   void RandomVectorState::initNLKernel(const bool force)
   {
      this->mRandom.setResolution(this->spRes());
      IVectorEquation::initNLKernel(force);
   }

   void RandomVectorState::setIdentity(const std::size_t name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void RandomVectorState::setSpectrum(const FieldComponents::Spectral::Id comp, const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const SpecialId special)
   {
      this->mRandom.setSpectrum(comp, min, max, ratio1D, ratio2D, special);
   }

   void RandomVectorState::setSpectrum(const FieldComponents::Spectral::Id comp, const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const MHDFloat ratio3D, const SpecialId special)
   {
      this->mRandom.setSpectrum(comp, min, max, ratio1D, ratio2D, ratio3D, special);
   }

   void RandomVectorState::setCoupling()
   {
      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Source) = true;
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

   MHDVariant RandomVectorState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
   {
      return this->mRandom.sourceTerm(compId, i, j, k);
   }

   void RandomVectorState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::After::id());

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add unknown to requirements: is scalar?
      auto& req = this->mRequirements.addField(this->name(), FieldRequirement(false, ss.spectral(), ss.physical()));
      req.enableSpectral();
      req.enablePhysical();
   }

   void RandomVectorState::setNLComponents()
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
