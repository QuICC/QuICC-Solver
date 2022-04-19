/**
 * @file TiltedScalarFieldVisualizer.cpp
 * @brief Source of the implementation of the tilted scalar field visualizer
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/Visualizers/TiltedScalarFieldVisualizer.hpp"

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/NonDimensional/Theta.hpp"
#include "QuICC/NonDimensional/Lower1D.hpp"
#include "QuICC/NonDimensional/Upper1D.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/TypeSelectors/TransformSelector.hpp"
#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Backward/P.hpp"

namespace QuICC {

namespace Equations {

   TiltedScalarFieldVisualizer::TiltedScalarFieldVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme)
      : IScalarEquation(spEqParams,spScheme), mViewField(true), mViewGradient(false)
   {
   }

   TiltedScalarFieldVisualizer::~TiltedScalarFieldVisualizer()
   {
   }

   void TiltedScalarFieldVisualizer::setIdentity(const std::size_t name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void TiltedScalarFieldVisualizer::setDataField(const std::size_t name)
   {
      this->mDataField = name;
   }

   void TiltedScalarFieldVisualizer::setFields(const bool viewField, const bool viewGradient)
   {
      this->mViewField = viewField;

      this->mViewGradient = viewGradient;
   }

   void TiltedScalarFieldVisualizer::setCoupling()
   {
      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Nonlinear) = true;

      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::WRAPPER, 0, features);
   }

   void TiltedScalarFieldVisualizer::computeNonlinear(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Get paramters
      MHDFloat eta3 = std::cos((Math::PI/180.)*this->eqParams().nd(NonDimensional::Theta::id()));
      MHDFloat eta2 = std::sin((Math::PI/180.)*this->eqParams().nd(NonDimensional::Theta::id()));

      // Initialize FFT
      Transform::Fft::MixedFftSelector   transform;
      Transform::Fft::MixedFftSelector::SharedSetupType spSetup = std::static_pointer_cast<Transform::Fft::MixedFftSelector::SetupType>(this->res().spTransformSetup(Dimensions::Transform::TRA3D));
      transform.init(spSetup);

      // Compute forward transform
      MatrixZ tmp(spSetup->bwdSize(), spSetup->blockSize());
      transform.forward(tmp, this->scalar(this->mDataField).dom(0).phys().data(), Transform::Forward::P::id());

      // Get Z grid
      int nK = this->res().sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nJ;
      int nI;
      MHDFloat zi = this->eqParams().nd(NonDimensional::Lower1D::id());
      MHDFloat zo = this->eqParams().nd(NonDimensional::Upper1D::id());
      Array gK = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nK, zi, zo);

      MHDFloat k_;
      int m = 0;
      nK = this->res().cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      for(int iK = 0; iK < nK; ++iK)
      {
         k_ = (1.0 - gK(this->res().cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iK)))/2.0;
         nJ = this->res().cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iK);
         for(int iJ = 0; iJ < nJ; ++iJ)
         {
            nI = this->res().cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DATB1D>(iK);
            for(int iI = 0; iI < nI; ++iI)
            {
               tmp(iI, m) = std::exp(MHDComplex(0.0, -k_*iI*(eta2/eta3)))*tmp(iI,m);
            }
            m++;
         }
      }

      // Compute backward transform
      transform.backward(rNLComp.rData(), tmp, Transform::Backward:P::id());
   }

   void TiltedScalarFieldVisualizer::useNonlinear(const Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id compId)
   {
      this->rUnknown().rDom(0).rPhys().rData() = rNLComp.data();
   }

   void TiltedScalarFieldVisualizer::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::After::id());

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add unknown to requirements: is scalar?
      auto& req = this->mRequirements.addField(this->name(), FieldRequirement(true, ss.spectral(), ss.physical()));
      req.enableSpectral();
      if(this->mViewField) req.enablePhysical();
      if(this->mViewGradient) req.enableGradient();

      // Add unknown to requirements: is scalar?
      auto& reqData = this->mRequirements.addField(this-mDataFiel, FieldRequirement(true, ss.spectral(), ss.physical()));
      reqData.enableSpectral();
      if(this->mViewField) reqData.enablePhysical();
      if(this->mViewGradient) reqData.enableGradient();
   }

}
}
