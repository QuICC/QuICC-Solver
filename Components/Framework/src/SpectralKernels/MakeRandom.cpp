/**
 * @file MakeRandom.cpp
 * @brief Source of trivial kernel to make field random with given scale
 */

// Configuration includes
//

// System includes
//
#include <time.h>

// External includes
//

// Class include
//
#include "QuICC/SpectralKernels/MakeRandom.hpp"

// Project includes
//
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Spectral {

namespace Kernel {

   MakeRandom::MakeRandom(const bool isComplex)
      : ISpectralKernel(isComplex), mStartSeed(1), mMin(-1.0), mMax(1.0)
   {
      timespec t;
      clock_gettime(CLOCK_REALTIME, &t);
      srand(t.tv_sec + t.tv_nsec);

      this->mStartSeed = t.tv_sec + t.tv_nsec;
   }

   MakeRandom::~MakeRandom()
   {
   }

   void MakeRandom::init(const MHDFloat minVal, const MHDFloat maxVal)
   {
      if(maxVal <= minVal)
      {
         throw std::logic_error("Incompatible random spectrum properties requested!");
      }

      this->mMin = minVal;
      this->mMax = maxVal;

      // Default spectrum ratio to 1.0
      if(this->mRatio.empty())
      {
         this->mRatio.push_back(1.0);
         this->mRatio.push_back(1.0);
         this->mRatio.push_back(1.0);
      }

      // Default number of zero modes at end of spectrum
      if(this->mZeros.empty())
      {
         this->mZeros.push_back(4);
         this->mZeros.push_back(4);
         this->mZeros.push_back(4);
      }
   }

   void MakeRandom::setRatio(const std::vector<MHDFloat>& ratios)
   {
      if(*std::min_element(ratios.begin(), ratios.end()) < 1.0)
      {
         throw std::logic_error("Incompatible random spectrum properties requested!");
      }

      this->mRatio = ratios;
   }

   void MakeRandom::setZeros(const std::vector<int>& zeros)
   {
      if(*std::min_element(zeros.begin(), zeros.end()) < 0)
      {
         throw std::logic_error("Number of zeros at end of spectrum should be positive!");
      }

      this->mZeros = zeros;
   }

   MHDFloat MakeRandom::scalingRatio(const int i, const int n, const int dim) const
   {
      if(static_cast<std::size_t>(dim) > this->mRatio.size())
      {
         throw std::logic_error("Scaling ratio is out of dimension bounds");
      }

      return std::exp(-static_cast<MHDFloat>(i)*std::log(this->mRatio.at(dim))/static_cast<MHDFloat>(n));
   }

   MHDVariant MakeRandom::compute(const int i, const int j, const int k) const
   {
      if(this->mIsComplex)
      {
         MHDComplex tmp;
         this->mode(tmp, i, j, k);
         return tmp;
      } else
      {
         MHDFloat tmp;
         this->mode(tmp, i, j, k);
         return tmp;
      }
   }

   void MakeRandom::value(MHDFloat& val, const unsigned int seed) const
   {
      // Change seed
      if(seed != 1)
      {
         srand(this->mStartSeed + seed);
      }

      // Get random value in range
      val = ((this->mMin-this->mMax)*static_cast<MHDFloat>(rand())/RAND_MAX)+this->mMax;

      // Change seed
      if(seed != 1)
      {
         timespec t;
         clock_gettime(CLOCK_REALTIME, &t);
         srand(t.tv_sec + t.tv_nsec);
      }
   }

   void MakeRandom::value(MHDComplex& val, const unsigned int seed) const
   {
      MHDFloat tmp;
      this->value(tmp, seed);
      val.real(tmp);
      this->value(tmp, seed);
      val.imag(tmp);
   }

   template <typename T> void MakeRandom::mode(T& val, const int i, const int j, const int k) const
   {
      // Initalize to zero
      val = T(0.0);

      if(this->res().sim().ss().dimension() == 3)
      {
         // Get first dimension
         const int n1D = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
         // Get second dimension
         int n2D = this->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);
         // Get third dimension
         const int n3D = this->res().sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

         // Get simulation wide indexes
         int j_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
         int k_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

         const int z1D = this->mZeros.at(0);
         int z2D = this->mZeros.at(1);
         const int z3D = this->mZeros.at(2);
         if(this->res().sim().ss().has(SpatialScheme::Feature::FourierIndex23))
         {
//            if(this->mSpecial == ONLYMEAN && !(j_ == 0 && k_ == 0))
//            {
//               val = T(0);
//            } else if(this->mSpecial == NOMEAN && j_ == 0 && k_ == 0)
//            {
//               val = T(0);
//            }

            z2D = z2D/2;
            if(k_ >= n2D/2)
            {
               k_ = n2D - k_;
            }
            n2D = n2D/2;
         }

         if(i < n1D-z1D && j_ < n3D - z3D && k_ < n2D - z2D)
         {
            // Compute scaling factors
            MHDFloat a1D = this->scalingRatio(i, n1D, 0);
            MHDFloat a2D = this->scalingRatio(k_, n2D, 1);
            MHDFloat a3D = this->scalingRatio(j_, n3D, 2);

            // Generate random value
            this->value(val);

            // Scale value
            val *= a1D*a2D*a3D;
         }
      } else if(this->res().sim().ss().dimension() == 2)
      {
         // Get first dimension
         const int n1D = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
         // Get second dimension
         const int n2D = this->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

         // Get simulation wide indexes
         int j_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j);

         const int z1D = this->mZeros.at(0);
         const int z2D = this->mZeros.at(1);

         if(i < n1D-z1D && j_ < n2D - z2D)
         {
            // Compute scaling factors
            MHDFloat a1D = this->scalingRatio(i, n1D, 0);
            MHDFloat a2D = this->scalingRatio(j_, n2D, 1);

            // Generate random value
            this->value(val);

            // Scale value
            val *= a1D*a2D;
         }
      } else
      {
         // Get first dimension
         const int n1D = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

         const int z1D = this->mZeros.at(0);

         if(i < n1D-z1D)
         {
            // Compute scaling factors
            MHDFloat a1D = this->scalingRatio(i, n1D, 0);

            // Generate random value
            this->value(val);

            // Scale value
            val *= a1D;
         }
      }

      this->constrain(val, i, j, k);
   }

   void MakeRandom::constrain(MHDFloat& val, const int i, const int j, const int k) const
   {
      // Do nothing for the moment
   }

   void MakeRandom::constrain(MHDComplex& val, const int i, const int j, const int k) const
   {
      // Force real value for k = 0 Fourier modes
      if(this->res().sim().ss().has(SpatialScheme::Feature::FourierIndex3))
      {
         if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 0)
         {
            val.imag(0.0);
         }
      // Force real value for k = 0
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

            MHDFloat tmp;
            this->value(tmp, seed);
            val.real(tmp);
            this->value(tmp, seed+3);
            if(k2D < n2D/2)
            {
               val.imag(tmp);
            } else
            {
               val.imag(-tmp);
            }

         } else if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 0)
         {
            val.imag(0.0);
         }
      } else if(this->res().sim().ss().has(SpatialScheme::Feature::FourierIndex2))
      {
         if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k) == 0)
         {
            val.imag(0.0);
         }
      }

      // Assume spherical harmonic for shell and shere
      if(this->res().sim().ss().has(SpatialScheme::Feature::SphereGeometry) || this->res().sim().ss().has(SpatialScheme::Feature::SphereGeometry))
      {
         if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 0)
         {
            val = 0.0;
         }
      }
   }

}
}
}
