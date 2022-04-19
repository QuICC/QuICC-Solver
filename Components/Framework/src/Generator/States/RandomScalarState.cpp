/**
 * @file RandomScalarState.cpp
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
#include "QuICC/Generator/States/RandomScalarState.hpp"

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/SolveTiming/Before.hpp"

namespace QuICC {

namespace Equations {

   RandomScalarState::RandomScalarState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IScalarEquation(spEqParams,spScheme,spBackend)
   {
   }

   RandomScalarState::~RandomScalarState()
   {
   }

   void RandomScalarState::initNLKernel(const bool force)
   {
      this->mRandom.setResolution(this->spRes());
      IScalarEquation::initNLKernel(force);
   }

   void RandomScalarState::setIdentity(const std::size_t name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void RandomScalarState::setSpectrum(const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const SpecialId special)
   {
      this->mRandom.setSpectrum(FieldComponents::Spectral::SCALAR, min, max, ratio1D, ratio2D, special);
   }

   void RandomScalarState::setSpectrum(const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const MHDFloat ratio3D, const SpecialId special)
   {
      this->mRandom.setSpectrum(FieldComponents::Spectral::SCALAR, min, max, ratio1D, ratio2D, ratio3D, special);
   }

   void RandomScalarState::setCoupling()
   {
      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Source) = true;
      features.at(CouplingFeature::AllowExplicit) = false;

      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, features);
   }

   MHDVariant RandomScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      return this->mRandom.sourceTerm(compId, i, j, k);
   }

   void RandomScalarState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::Before::id());

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add unknown to requirements: is scalar?
      auto& req = this->mRequirements.addField(this->name(), FieldRequirement(true, ss.spectral(), ss.physical()));
      req.enableSpectral();
      req.enablePhysical();
   }

}
}
