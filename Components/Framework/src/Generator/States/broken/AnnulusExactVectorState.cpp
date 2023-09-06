/**
 * @file AnnulusExactVectorState.cpp
 * @brief Source of the implementation of the equation to generate an exact vector solution in an annulus
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
#include "QuICC/Generator/States/AnnulusExactVectorState.hpp"

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "QuICC/NonDimensional/Lower3d.hpp"
#include "QuICC/NonDimensional/Upper3d.hpp"

namespace QuICC {

namespace Equations {

   AnnulusExactVectorState::AnnulusExactVectorState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme)
      : IVectorEquation(spEqParams,spScheme)
   {
   }

   AnnulusExactVectorState::~AnnulusExactVectorState()
   {
   }

   void AnnulusExactVectorState::setIdentity(const std::size_t name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void AnnulusExactVectorState::setStateType(const FieldComponents::Physical::Id compId, const AnnulusExactStateIds::Id id)
   {
      this->mTypeId.insert(std::make_pair(compId, id));
   }

   void AnnulusExactVectorState::setModeOptions(const FieldComponents::Physical::Id compId, const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2, const MHDFloat a3, const MHDFloat k3)
   {
      Array tmp(3);
      tmp(0) = a1;
      tmp(1) = a2;
      tmp(2) = a3;
      this->mModeA.insert(std::make_pair(compId, tmp));

      tmp(0) = k1;
      tmp(1) = k2;
      tmp(2) = k3;
      this->mModeK.insert(std::make_pair(compId, tmp));
   }

   void AnnulusExactVectorState::setCoupling()
   {
      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Nonlinear) = true;
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

   void AnnulusExactVectorState::computeNonlinear(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id compId) const
   {
      AnnulusExactStateIds::Id typeId = this->mTypeId.find(compId)->second;
      Array modeA = this->mModeA.find(compId)->second;
      Array modeK = this->mModeK.find(compId)->second;

      if(typeId == AnnulusExactStateIds::CONSTANT)
      {
         rNLComp.rData().setConstant(modeA(0)*modeA(1)*modeA(2));
      } else
      {
         int nR = this->res().sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nT = this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nZ = this->res().sim().dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array& gR = this->mspMesh->at(0);
         Array& gT = this->mspMesh->at(1);
         Array& gZ = this->mspMesh->at(2);

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

                  if(typeId == AnnulusExactStateIds::POLYCOSPOLY)
                  {
                     valR = AnnulusExactStateIds::poly(modeA(0),modeK(0),r_);
                     valT = AnnulusExactStateIds::cos(modeA(1),modeK(1),t_);
                     valZ = AnnulusExactStateIds::poly(modeA(2),modeK(2),z_);

                  } else if(typeId == AnnulusExactStateIds::POLYSINPOLY)
                  {
                     valR = AnnulusExactStateIds::poly(modeA(0),modeK(0),r_);
                     valT = AnnulusExactStateIds::sin(modeA(1),modeK(1),t_);
                     valZ = AnnulusExactStateIds::poly(modeA(2),modeK(2),z_);

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

   MHDVariant AnnulusExactVectorState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
   {
      return MHDComplex(0);
   }

   void AnnulusExactVectorState::setRequirements()
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

   void AnnulusExactVectorState::setNLComponents()
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
