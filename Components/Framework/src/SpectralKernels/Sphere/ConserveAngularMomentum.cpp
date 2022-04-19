/**
 * @file ConserveAngularMomentum.cpp
 * @brief Source of angular momentum conservation kernel in a sphere
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SpectralKernels/Sphere/ConserveAngularMomentum.hpp"

// Project includes
//
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Polynomial/Worland/Operators.hpp"

namespace QuICC {

namespace Spectral {

namespace Kernel {

namespace Sphere {

   ConserveAngularMomentum::ConserveAngularMomentum(const bool isComplex)
      : ISpectralKernel(isComplex), mHasM0(false), mHasM1(false), mM0j(-1), mM0k(-1), mM1j(-1), mM1k(-1)
   {
   }

   ConserveAngularMomentum::~ConserveAngularMomentum()
   {
   }

   void ConserveAngularMomentum::init(const bool hasMOrdering)
   {
      if(hasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int m_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

               if(l_ == 1 && m_ == 0)
               {
                  this->mM0j = j;
                  this->mM0k = k;
                  this->mHasM0 = true;
               } else if(l_ == 1 && m_ == 1)
               {
                  this->mM1j = j;
                  this->mM1k = k;
                  this->mHasM1 = true;
               }
            }
         }
      } else
      {
         // Loop over harmonic degree l
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int m_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
               if(l_ == 1 && m_ == 0)
               {
                  this->mM0j = j;
                  this->mM0k = k;
                  this->mHasM0 = true;
               } else if(l_ == 1 && m_ == 1)
               {
                  this->mM1j = j;
                  this->mM1k = k;
                  this->mHasM1 = true;
               }
            }
         }
      }

      // Compute operator if required
      if(this->mHasM0 || this->mHasM1)
      {
         int nN = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
         internal::Matrix iop;
         Polynomial::Worland::Operators::integrateRpWnl(iop, 1, 3, nN);
         iop.col(0).bottomRows(iop.rows()-1) /= -iop(0,0);
         iop(0,0) = MHD_MP(0);
         this->mOp = iop.cast<MHDFloat>();
         assert(this->mOp.rows() == nN && this->mOp.cols() == 1);
      }
   }

   MHDVariant ConserveAngularMomentum::compute(const int i, const int j, const int k) const
   {
      if(this->mIsComplex)
      {
         return MHDComplex(0,0);
      } else
      {
         return 0.0;
      }
   }

   void ConserveAngularMomentum::apply()
   {
      assert(this->mScalars.size() == 0);
      assert(this->mVectors.size() == 1);
      auto& field = this->mVectors.begin()->second;

      ArrayZ mom;
      if(this->mHasM1)
      {
         std::visit([&](auto&& f){
            mom = (this->mOp.transpose()*f->dom(0).total().comp(FieldComponents::Spectral::TOR).profile(this->mM1j, this->mM1k));
            f->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).setPoint(mom(0), 0, this->mM1j,this->mM1k);
         }, field);
      }

      if(this->mHasM0)
      {
         std::visit([&](auto&& f){
            mom = (this->mOp.transpose()*f->dom(0).total().comp(FieldComponents::Spectral::TOR).profile(this->mM0j, this->mM0k));
            f->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).setPoint(mom(0), 0, this->mM0j, this->mM0k);
         }, field);
      }
   }

}
}
}
}
