/**
 * @file AnnulusExactScalarState.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution in an annulus
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
#include "QuICC/Generator/States/AnnulusExactScalarState.hpp"

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/NonDimensional/Lower1D.hpp"
#include "QuICC/NonDimensional/Upper1D.hpp"
#include "QuICC/NonDimensional/Lower3D.hpp"
#include "QuICC/NonDimensional/Upper3D.hpp"

namespace QuICC {

namespace Equations {

   AnnulusExactScalarState::AnnulusExactScalarState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme)
      : IScalarEquation(spEqParams,spScheme), mTypeId(AnnulusExactStateIds::CONSTANT), mModeA(3), mModeK(3)
   {
   }

   AnnulusExactScalarState::~AnnulusExactScalarState()
   {
   }

   void AnnulusExactScalarState::setIdentity(const std::size_t name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void AnnulusExactScalarState::setStateType(const AnnulusExactStateIds::Id id)
   {
      this->mTypeId = id;
   }

   void AnnulusExactScalarState::setModeOptions(const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2, const MHDFloat a3, const MHDFloat k3)
   {
      this->mModeA(0) = a1;
      this->mModeA(1) = a2;
      this->mModeA(2) = a3;

      this->mModeK(0) = k1;
      this->mModeK(1) = k2;
      this->mModeK(2) = k3;
   }

   void AnnulusExactScalarState::setCoupling()
   {
      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Nonlinear) = true;
      features.at(CouplingFeature::AllowExplicit) = false;

      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, features);
   }

   void AnnulusExactScalarState::computeNonlinear(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Physical::SCALAR);

      if(this->mTypeId == AnnulusExactStateIds::CONSTANT)
      {
         rNLComp.rData().setConstant(this->mModeA(0)*this->mModeA(1)*this->mModeA(2));
      } else
      {
         int nR = this->res().sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nT = this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nZ = this->res().sim().dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array& gR = this->mspMesh->at(0);
         Array& gT = this->mspMesh->at(1);
         Array gZ = this->mspMesh->at(2);

         MHDFloat r_;
         MHDFloat t_;
         MHDFloat z_;

         nR = this->res().cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         for(int iR = 0; iR < nR; ++iR)
         {
            r_ = gR(this->res().cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            nT = this->res().cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR);
            for(int iT = 0; iT < nT; ++iT)
            {
               t_ = gT(this->res().cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iT, iR));
               for(int iZ = 0; iZ < nZ; ++iZ)
               {
                  z_ = gZ(iZ);

                  MHDFloat valZ = 0.0;
                  MHDFloat valT = 0.0;
                  MHDFloat valR = 0.0;

                  if(this->mTypeId == AnnulusExactStateIds::POLYCOSPOLY)
                  {
                     valR = AnnulusExactStateIds::poly(this->mModeA(0),this->mModeK(0),r_);
                     valT = AnnulusExactStateIds::cos(this->mModeA(1),this->mModeK(1),t_);
                     valZ = AnnulusExactStateIds::poly(this->mModeA(2),this->mModeK(2),z_);

                  } else if(this->mTypeId == AnnulusExactStateIds::POLYSINPOLY)
                  {
                     valR = AnnulusExactStateIds::poly(this->mModeA(0),this->mModeK(0),r_);
                     valT = AnnulusExactStateIds::sin(this->mModeA(1),this->mModeK(1),t_);
                     valZ = AnnulusExactStateIds::poly(this->mModeA(2),this->mModeK(2),z_);

                  } else
                  {
                     throw std::logic_error("Unknown exact state");
                  }

                  rNLComp.setPoint(valR*valT*valZ, iZ, iT, iR);
               }
            }
         }
      }
   }

    MHDVariant AnnulusExactScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      return MHDComplex(0);
    }

   void AnnulusExactScalarState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::After::id());

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add unknown to requirements: is scalar?
      auto& req = this->mRequirements.addField(this->name(), FieldRequirement(true, ss.spectral(), ss.physical()));
      req.enableSpectral();
   }

}
}
