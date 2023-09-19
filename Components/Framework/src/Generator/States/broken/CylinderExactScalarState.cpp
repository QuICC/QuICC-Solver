/**
 * @file CylinderExactScalarState.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution in a cylinder
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
#include "QuICC/Generator/States/CylinderExactScalarState.hpp"

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Constants.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/NonDimensional/Lower3d.hpp"
#include "QuICC/NonDimensional/Upper3d.hpp"

namespace QuICC {

namespace Equations {

   CylinderExactScalarState::CylinderExactScalarState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme)
      : IScalarEquation(spEqParams,spScheme), mTypeId(CylinderExactStateIds::NOTUSED), mSpecTypeId(CylinderExactStateIds::NOTUSED), mModeA(3), mModeK(3)
   {
   }

   CylinderExactScalarState::~CylinderExactScalarState()
   {
   }

   void CylinderExactScalarState::setIdentity(const std::size_t name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void CylinderExactScalarState::setPhysicalType(const CylinderExactStateIds::Id id)
   {
      this->mTypeId = id;
   }

   void CylinderExactScalarState::setSpectralType(const CylinderExactStateIds::Id id)
   {
      this->mSpecTypeId = id;
   }

   void CylinderExactScalarState::setModeOptions(const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2, const MHDFloat a3, const MHDFloat k3)
   {
      this->mModeA(0) = a1;
      this->mModeA(1) = a2;
      this->mModeA(2) = a3;

      this->mModeK(0) = k1;
      this->mModeK(1) = k2;
      this->mModeK(2) = k3;
   }

   void CylinderExactScalarState::setCoupling()
   {
      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Nonlinear) = (this->mTypeId != CylinderExactStateIds::NOTUSED);
      features.at(CouplingFeature::Source) = (this->mSpecTypeId != CylinderExactStateIds::NOTUSED);
      features.at(CouplingFeature::AllowExplicit) = false;

      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, features);
   }

   void CylinderExactScalarState::computeNonlinear(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Physical::SCALAR);

      if(this->mTypeId == CylinderExactStateIds::CONSTANT)
      {
         rNLComp.rData().setConstant(this->mModeA(0)*this->mModeA(1)*this->mModeA(2));
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

                  if(this->mTypeId == CylinderExactStateIds::POLYCOSPOLY)
                  {
                     MHDFloat valR = CylinderExactStateIds::poly(this->mModeA(0),this->mModeK(0),r_);
                     MHDFloat valT = CylinderExactStateIds::cos(this->mModeA(1),this->mModeK(1),t_);
                     MHDFloat valZ = CylinderExactStateIds::poly(this->mModeA(2),this->mModeK(2),z_);

                     rNLComp.setPoint(valR*valT*valZ, iZ, iT, iR);

                  } else if(this->mTypeId == CylinderExactStateIds::POLYSINPOLY)
                  {
                     MHDFloat valR = CylinderExactStateIds::poly(this->mModeA(0),this->mModeK(0),r_);
                     MHDFloat valT = CylinderExactStateIds::sin(this->mModeA(1),this->mModeK(1),t_);
                     MHDFloat valZ = CylinderExactStateIds::poly(this->mModeA(2),this->mModeK(2),z_);

                     rNLComp.setPoint(valR*valT*valZ, iZ, iT, iR);

                  } else
                  {
                     throw std::logic_error("Unknown exact state");
                  }
               }
            }
         }
      }
   }

    MHDVariant CylinderExactScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int iR, const int iZ, const int iM) const
    {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      if(this->mSpecTypeId == CylinderExactStateIds::SPEC_UNIT)
      {
         int m = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(iM);
         if(m == 1)
         {
            int z = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(iZ,iM);
            if(z == 0)
            {
               if(iR < 6)
               {
                  return MHDComplex(1.0);
               }
            }
         }

         return MHDComplex(0);
      } else
      {
         return MHDComplex(0);
      }
    }

   void CylinderExactScalarState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::After::id());

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add unknown to requirements: is scalar?
      auto& req = this->mRequirements.addField(this->name(), FieldRequirement(true, ss.spectral(), ss.physical()));
      req.enableSpectral();
      req.enablePhysical();
   }

}
}
