/**
 * @file TorPolHarmonic.cpp
 * @brief Source of IMP Magnetic state from Matsui generator kernel
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/States/Kernels/Sphere/TorPolHarmonic.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"
#include "QuICC/SpectralKernels/Typedefs.hpp"
#include "QuICC/Generator/States/Kernels/Tools/SphericalHarmonic.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Sphere {

   TorPolHarmonic::TorPolHarmonic()
      : IPhysicalKernel()
   {
   }

   TorPolHarmonic::~TorPolHarmonic()
   {
   }

   void TorPolHarmonic::setModes(const FieldComponents::Spectral::Id compId, const Spectral::Kernel::Complex3DMapType& modes)
   {
      auto spModes = std::make_shared<Spectral::Kernel::Complex3DMapType>();
      *spModes = modes;

      this->mSHModes.insert(std::make_pair(compId, spModes));
   }

   void TorPolHarmonic::init()
   {
   }

   void TorPolHarmonic::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      Profiler::RegionFixture<2> fix("TorPolHarmonicCompute");

      // Initialize to zero
      rNLComp.rData().setZero();

      int nR = this->spRes()->sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nTh = this->spRes()->sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nPh = this->spRes()->sim().dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      Array& rGrid = this->mspMesh->at(0);
      Array& thGrid = this->mspMesh->at(1);
      Array& phGrid = this->mspMesh->at(2);

      typedef Spectral::Kernel::Complex3DMapType::const_iterator ModeIt;
      std::pair<ModeIt, ModeIt>  modeRange;

      if(!this->mSHModes.find(FieldComponents::Spectral::TOR)->second->empty())
      {
         Array funcSH(nPh);
         modeRange.first = this->mSHModes.find(FieldComponents::Spectral::TOR)->second->begin();
         modeRange.second = this->mSHModes.find(FieldComponents::Spectral::TOR)->second->end();

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

               for(ModeIt it = modeRange.first; it != modeRange.second; ++it)
               {
                  int l = it->first.first;
                  int m = it->first.second;

                  MHDComplex ct = 0.0;
                  for(auto rIt = it->second.begin(); rIt != it->second.end(); ++rIt)
                  {
                     MHDComplex tmp;

                     int k = l + 2*rIt->first;
                     tmp.real(rIt->second.real()*std::pow(r,k));
                     tmp.imag(rIt->second.imag()*std::pow(r,k));

                     ct += tmp;
                  }

                  Tools::SphericalHarmonic::Torlm(funcSH, l, m, id, ct, theta, phGrid);

                  rNLComp.addProfile(funcSH,iTh,iR);
               }
            }
         }
      }

      if(!this->mSHModes.find(FieldComponents::Spectral::POL)->second->empty())
      {
         Array funcSH(nPh);

         modeRange.first = this->mSHModes.find(FieldComponents::Spectral::POL)->second->begin();
         modeRange.second = this->mSHModes.find(FieldComponents::Spectral::POL)->second->end();

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
               for(ModeIt it = modeRange.first; it != modeRange.second; ++it)
               {
                  int l = it->first.first;
                  int m = it->first.second;

                  MHDComplex cq = 0.0;
                  MHDComplex cs = 0.0;
                  for(auto rIt = it->second.begin(); rIt != it->second.end(); ++rIt)
                  {
                     MHDComplex tmpq, tmps;

                     int k = l + 2*rIt->first - 1;
                     tmpq.real(rIt->second.real()*std::pow(r,k));
                     tmpq.imag(rIt->second.imag()*std::pow(r,k));

                     MHDFloat dk = static_cast<MHDFloat>(k+2);
                     tmps.real(rIt->second.real()*dk*std::pow(r,k));
                     tmps.imag(rIt->second.imag()*dk*std::pow(r,k));

                     cq += tmpq;
                     cs += tmps;
                  }

                  Tools::SphericalHarmonic::Pollm(funcSH, l, m, id, cq, cs, theta, phGrid);

                  rNLComp.addProfile(funcSH,iTh,iR);
               }
            }
         }
      }
   }

}
}
}
}
