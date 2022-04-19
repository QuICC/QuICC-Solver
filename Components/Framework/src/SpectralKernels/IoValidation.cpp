/**
 * @file IoValidation.cpp
 * @brief Source of trivial spectral kernel to generate IO validation field
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SpectralKernels/IoValidation.hpp"

// Project includes
//
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Spectral {

namespace Kernel {

   IoValidation::IoValidation(const bool isComplex)
      : ISpectralKernel(isComplex), mRange(100)
   {
   }

   IoValidation::~IoValidation()
   {
   }

   void IoValidation::init(const MHDFloat range)
   {
      this->mRange = range;
   }

   MHDVariant IoValidation::compute(const int i, const int j, const int k) const
   {
      // Create key
      std::pair<int,int>   key;
      if(this->spRes()->sim().ss().has(SpatialScheme::Feature::SpectralOrdering132))
      {
         key = std::make_pair(this->spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k),this->spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k));
      } else
      {
         key = std::make_pair(this->spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k),this->spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k));
      }

      MHDFloat range = 1000;
      MHDFloat idVal = key.first*range*range + key.second*range + i;

      // Set value
      if(this->mIsComplex)
      {
         return MHDComplex(idVal, -idVal);
      } else
      {
         return idVal;
      }
   }

}
}
}
