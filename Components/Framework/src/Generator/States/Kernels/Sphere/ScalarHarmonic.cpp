/**
 * @file ScalarHarmonic.cpp
 * @brief Source of scalar harmonic state generator kernel
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/States/Kernels/Sphere/ScalarHarmonic.hpp"

// Project includes
//
#include "QuICC/Generator/States/Kernels/Tools/SphericalHarmonic.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Sphere {

   ScalarHarmonic::ScalarHarmonic()
      : IPhysicalKernel()
   {
   }

   ScalarHarmonic::~ScalarHarmonic()
   {
   }

   void ScalarHarmonic::init(const Spectral::Kernel::Complex3DMapType& shModes)
   {
      this->mSHModes = std::make_shared<Spectral::Kernel::Complex3DMapType>();
      *this->mSHModes = shModes;
   }

   void ScalarHarmonic::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Initialize to zero
      rNLComp.rData().setZero();

      int nR = this->spRes()->sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nTh = this->spRes()->sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nPh = this->spRes()->sim().dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      Array& rGrid = this->mspMesh->at(0);
      Array& thGrid = this->mspMesh->at(1);
      Array& phGrid = this->mspMesh->at(2);

      Array funcSH(nPh);
      typedef Spectral::Kernel::Complex3DMapType::const_iterator ModeIt;
      //typedef Spectral::Kernel::ComplexFastMapType::const_iterator RadialIt;
      ModeIt it;
      std::pair<ModeIt, ModeIt>  modeRange = std::make_pair(this->mSHModes->begin(), this->mSHModes->end());

      MHDFloat r;
      MHDFloat theta;
      nR = this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      for(int iR = 0; iR < nR; ++iR)
      {
         r = rGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
         nTh = this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR);
         for(int iTh = 0; iTh < nTh; ++iTh)
         {
            theta = thGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));
            for(it = modeRange.first; it != modeRange.second; ++it)
            {
               int l = it->first.first;
               int m = it->first.second;

               MHDComplex c = 0.0;
               for(auto rIt = it->second.begin(); rIt != it->second.end(); ++rIt)
               {
                  MHDComplex tmp;

                  int k = l + 2*rIt->first;
                  tmp.real(rIt->second.real()*std::pow(r,k));
                  tmp.imag(rIt->second.imag()*std::pow(r,k));

                  c += tmp;
               }

               Tools::SphericalHarmonic::Ylm(funcSH, l, m, c, theta, phGrid);

               rNLComp.addProfile(funcSH,iTh,iR);
            }
         }
      }
   }

}
}
}
}
