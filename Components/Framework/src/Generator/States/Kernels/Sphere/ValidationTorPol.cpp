/**
 * @file ValidationTorPol.cpp
 * @brief Source of benchmark state C2 for velocity field generator kernel
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/States/Kernels/Sphere/ValidationTorPol.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Sphere {

   ValidationTorPol::ValidationTorPol()
      : IPhysicalKernel(), mC(1.0)
   {
   }

   ValidationTorPol::~ValidationTorPol()
   {
   }

   void ValidationTorPol::init(const MHDFloat c)
   {
      this->mC = c;
   }

   Array ValidationTorPol::x(const int iR, const int iTh) const
   {
      Array& rGrid = this->mspMesh->at(0);
      Array& thGrid = this->mspMesh->at(1);
      Array& phGrid = this->mspMesh->at(2);
      MHDFloat r = rGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
      MHDFloat theta = thGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

      return r*std::sin(theta)*(phGrid).array().cos();
   }

   Array ValidationTorPol::y(const int iR, const int iTh) const
   {
      Array& rGrid = this->mspMesh->at(0);
      Array& thGrid = this->mspMesh->at(1);
      Array& phGrid = this->mspMesh->at(2);
      MHDFloat r = rGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
      MHDFloat theta = thGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

      return r*std::sin(theta)*(phGrid).array().sin();
   }

   Array ValidationTorPol::z(const int iR, const int iTh) const
   {
      Array& rGrid = this->mspMesh->at(0);
      Array& thGrid = this->mspMesh->at(1);
      MHDFloat r = rGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
      MHDFloat theta = thGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

      return r*std::cos(theta)*Array::Ones(this->spRes()->sim().dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL));
   }

   Array ValidationTorPol::vx(const int iR, const int iTh) const
   {
      return this->mC*this->x(iR, iTh).array().pow(2)*this->y(iR, iTh).array().pow(2)*this->z(iR, iTh).array().pow(2);
   }

   Array ValidationTorPol::vy(const int iR, const int iTh) const
   {
      return -this->mC*this->x(iR, iTh).array()*this->y(iR, iTh).array().pow(3)*this->z(iR, iTh).array().pow(2);
   }

   Array ValidationTorPol::vz(const int iR, const int iTh) const
   {
      return this->mC*(1./3.)*this->x(iR, iTh).array()*this->y(iR, iTh).array().pow(2)*this->z(iR, iTh).array().pow(3);
   }

   void ValidationTorPol::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Initialize to zero
      rNLComp.rData().setZero();

      int nR = this->spRes()->sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nTh = this->spRes()->sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      //int nPh = this->spRes()->sim().dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      //Array& rGrid = this->mspMesh->at(0);
      Array& thGrid = this->mspMesh->at(1);
      Array& phGrid = this->mspMesh->at(2);

      Array cosPhi = (phGrid).array().cos();
      Array sinPhi = (phGrid).array().sin();

      MHDFloat theta;
      nR = this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      for(int iR = 0; iR < nR; ++iR)
      {
         //MHDFloat r = rGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
         nTh = this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR);
         for(int iTh = 0; iTh < nTh; ++iTh)
         {
            theta = thGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

            if(id == FieldComponents::Physical::R)
            {
               Array vr = std::sin(theta)*cosPhi.array()*vx(iR, iTh).array() + std::sin(theta)*sinPhi.array()*vy(iR, iTh).array() + std::cos(theta)*vz(iR, iTh).array();
               rNLComp.addProfile(vr,iTh,iR);
            } else if(id == FieldComponents::Physical::THETA)
            {
               Array vtheta = std::cos(theta)*cosPhi.array()*vx(iR, iTh).array() + std::cos(theta)*sinPhi.array()*vy(iR, iTh).array() - std::sin(theta)*vz(iR, iTh).array();
               rNLComp.addProfile(vtheta,iTh,iR);
            } else if(id == FieldComponents::Physical::PHI)
            {
               Array vphi = -sinPhi.array()*vx(iR, iTh).array() + cosPhi.array()*vy(iR, iTh).array();
               rNLComp.addProfile(vphi,iTh,iR);
            }
         }
      }

   }

}
}
}
}
