/**
 * @file ConserveAngularMomentum.cpp
 * @brief Source of angular momentum conservation kernel in a sphere
 */

// System includes
//

// Class include
//
#include "QuICC/SpectralKernels/Sphere/ConserveAngularMomentum.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Polynomial/Worland/Operators.hpp"
#include "QuICC/SolveTiming/After.hpp"

namespace QuICC {

namespace Spectral {

namespace Kernel {

namespace Sphere {

   ConserveAngularMomentum::ConserveAngularMomentum(const bool isComplex)
      : ISpectralKernel(isComplex), mHasM0(false), mHasM1(false), mM0j(-1), mM0k(-1), mM1j(-1), mM1k(-1)
   {
   }

   void ConserveAngularMomentum::init(const bool hasMOrdering)
   {
      const auto& tRes = *this->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
      if(hasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int m_ = tRes.idx<Dimensions::Data::DAT3D>(k);
            for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);

               if(l_ == 1 && m_ == 0)
               {
                  this->mM0j = j;
                  this->mM0k = k;
                  this->mHasM0 = true;
               }
               else if(l_ == 1 && m_ == 1)
               {
                  this->mM1j = j;
                  this->mM1k = k;
                  this->mHasM1 = true;
               }
            }
         }
      }
      else
      {
         // Loop over harmonic degree l
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l_ = tRes.idx<Dimensions::Data::DAT3D>(k);
            for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int m_ = tRes.idx<Dimensions::Data::DAT2D>(j,k);
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
      }
      else
      {
         return 0.0;
      }
   }

   void ConserveAngularMomentum::apply(const std::size_t timeId)
   {
      if(timeId == SolveTiming::After::id())
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

} // Sphere
} // Kernel
} // Spectral
} // QuICC
