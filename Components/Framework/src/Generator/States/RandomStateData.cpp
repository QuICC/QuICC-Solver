/**
 * @file RandomStateData.cpp
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
#include "QuICC/Generator/States/RandomStateData.hpp"

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Equations {

   RandomStateData::RandomStateData()
       : mStartSeed(1), mSpecial(NOTHING)
   {
      timespec t;
      clock_gettime(CLOCK_REALTIME, &t);
      srand(t.tv_sec + t.tv_nsec);

      this->mStartSeed = t.tv_sec + t.tv_nsec;
   }

   RandomStateData::~RandomStateData()
   {
   }

   void RandomStateData::setResolution(SharedResolution spRes)
   {
      this->mspRes = spRes;
   }

   const Resolution& RandomStateData::res() const
   {
      return *this->mspRes;
   }

   MHDFloat RandomStateData::scalingRatio(const FieldComponents::Spectral::Id comp, const int i, const int dim, const int n) const
   {
      MHDFloat ratio;
      switch(dim)
      {
         case 1:
            ratio = this->mRatio1D.find(comp)->second;
            break;
         case 2:
            ratio = this->mRatio2D.find(comp)->second;
            break;
         case 3:
            ratio = this->mRatio3D.find(comp)->second;
            break;
         default:
            throw std::logic_error("Unknown dimension requested");
            break;
      }

      return std::exp(-static_cast<MHDFloat>(i)*std::log(ratio)/static_cast<MHDFloat>(n));
   }

   void RandomStateData::setSpectrum(const FieldComponents::Spectral::Id comp, const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const SpecialId special)
   {
      if(max <= min || ratio1D < 1 || ratio2D < 1)
      {
         throw std::logic_error("Incompatible spectrum properties requested!");
      }

      // Set min and max range for values
      this->mMin.insert(std::make_pair(comp, min));
      this->mMax.insert(std::make_pair(comp, max));

      // Set spectrum ratios
      this->mRatio1D.insert(std::make_pair(comp, ratio1D));
      this->mRatio2D.insert(std::make_pair(comp, ratio2D));

      // Special ID
      this->mSpecial = special;
   }

   void RandomStateData::setSpectrum(const FieldComponents::Spectral::Id comp, const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const MHDFloat ratio3D, const SpecialId special)
   {
      if(max <= min || ratio1D < 1 || ratio2D < 1 || ratio3D < 1)
      {
         throw std::logic_error("Incompatible spectrum properties requested!");
      }

      // Set min and max range for values
      this->mMin.insert(std::make_pair(comp, min));
      this->mMax.insert(std::make_pair(comp, max));

      // Set spectrum ratios
      this->mRatio1D.insert(std::make_pair(comp, ratio1D));
      this->mRatio2D.insert(std::make_pair(comp, ratio2D));
      this->mRatio3D.insert(std::make_pair(comp, ratio3D));

      // Special ID
      this->mSpecial = special;
   }

   void RandomStateData::makeRandom(MHDFloat& val, const int i, const int j, const int k, const FieldComponents::Spectral::Id comp, const unsigned int seed) const
   {
      MHDFloat minVal = this->mMin.find(comp)->second;
      MHDFloat maxVal = this->mMax.find(comp)->second;

      if(seed == 1)
      {
         val = ((minVal-maxVal)*static_cast<MHDFloat>(rand())/RAND_MAX)+maxVal;
      } else
      {
         srand(this->mStartSeed + seed);
         val = ((minVal-maxVal)*static_cast<MHDFloat>(rand())/RAND_MAX)+maxVal;
         timespec t;
         clock_gettime(CLOCK_REALTIME, &t);
         srand(t.tv_sec + t.tv_nsec);
      }
   }

   void RandomStateData::makeRandom(MHDComplex& val, const int i, const int j, const int k, const FieldComponents::Spectral::Id comp) const
   {
      MHDFloat tmp;
      this->makeRandom(tmp, i, j, k, comp);

      val.real(tmp);

      if(this->res().sim().ss().has(SpatialScheme::Feature::FourierIndex3))
      {
         if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) != 0)
         {
            this->makeRandom(tmp, i, j, k, comp);
            val.imag(tmp);
         } else
         {
            val.imag(0.0);
         }
      } else if(this->res().sim().ss().has(SpatialScheme::Feature::FourierIndex23))
      {
         if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k) == 0 && this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) != 0)
         {
            unsigned int seed = 2;
            seed += this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATF1D>(i,k);

            int n2D = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
            int k2D = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            if(k2D < n2D/2)
            {
               seed += k2D;
            } else
            {
               seed += n2D - k2D;
            }

            this->makeRandom(tmp, i, j, k, comp, seed);
            val.real(tmp);
            this->makeRandom(tmp, i, j, k, comp, seed+3);
            if(k2D < n2D/2)
            {
               val.imag(tmp);
            } else
            {
               val.imag(-tmp);
            }

         } else if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) != 0)
         {
            this->makeRandom(tmp, i, j, k, comp);
            val.imag(tmp);
         } else
         {
            val.imag(0.0);
         }
      } else if(this->res().sim().ss().has(SpatialScheme::Feature::FourierIndex2))
      {
         if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k) != 0)
         {
            this->makeRandom(tmp, i, j, k, comp);
            val.imag(tmp);
         } else
         {
            val.imag(0.0);
         }
      }
   }

   MHDVariant RandomStateData::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
   {
      if(this->res().sim().ss().dimension() == 3)
      {
         // Get first dimension
         int n1D = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
         // Get second dimension
         int n2D = this->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);
         // Get third dimension
         int n3D = this->res().sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

         // Get simulation wide indexes
         int j_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
         int k_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

         int z1D = 4;
         int z2D = 4;
         int z3D = 4;
         if(this->res().sim().ss().has(SpatialScheme::Feature::FourierIndex23))
         {
            if(this->mSpecial == ONLYMEAN && !(j_ == 0 && k_ == 0))
            {
               return MHDComplex(0);
            } else if(this->mSpecial == NOMEAN && j_ == 0 && k_ == 0)
            {
               return MHDComplex(0);
            }

            z2D = 2;
            if(k_ >= n2D/2)
            {
               k_ = n2D - k_;
            }
            n2D = n2D/2;
         }

         if(i < n1D-z1D && j_ < n3D - z3D && k_ < n2D - z2D)
         {
            // Compute scaling factors
            MHDFloat a1D = this->scalingRatio(compId, i, 1, n1D);
            MHDFloat a2D = this->scalingRatio(compId, k_, 2, n2D);
            MHDFloat a3D = this->scalingRatio(compId, j_, 3, n3D);

            MHDComplex val;

            this->makeRandom(val, i, j, k, compId);

            return val*a1D*a2D*a3D;
         } else
         {
            return MHDComplex(0);
         }
      } else if(this->res().sim().ss().dimension() == 2)
      {
         // Get first dimension
         int n1D = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
         // Get second dimension
         int n2D = this->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

         // Get simulation wide indexes
         int j_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j);

         int z1D = 4;
         int z2D = 4;

         if(i < n1D-z1D && j_ < n2D - z2D)
         {
            // Compute scaling factors
            MHDFloat a1D = this->scalingRatio(compId, i, 1, n1D);
            MHDFloat a2D = this->scalingRatio(compId, j_, 2, n2D);

            MHDComplex val;

            this->makeRandom(val, i, j, k, compId);

            return val*a1D*a2D;
         } else
         {
            return MHDComplex(0);
         }
      } else
      {
         // Get first dimension
         int n1D = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

         int z1D = 4;

         if(i < n1D-z1D)
         {
            // Compute scaling factors
            MHDFloat a1D = this->scalingRatio(compId, i, 1, n1D);

            MHDComplex val;

            this->makeRandom(val, i, j, k, compId);

            return val*a1D;
         } else
         {
            return MHDComplex(0);
         }
      }
   }

}
}
