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
#include "QuICC/Enums/Dimensions.hpp"
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
      const auto& tRes = *this->spRes()->cpu()->dim(Dimensions::Transform::SPECTRAL);
      // Create key
      std::pair<int,int>   key;
      if(this->spRes()->sim().ss().has(SpatialScheme::Feature::SpectralOrdering132))
      {
         key = std::make_pair(tRes.idx<Dimensions::Data::DAT3D>(k),tRes.idx<Dimensions::Data::DAT2D>(j,k));
      } else
      {
         key = std::make_pair(tRes.idx<Dimensions::Data::DAT2D>(j,k),tRes.idx<Dimensions::Data::DAT3D>(k));
      }

      MHDFloat idVal = key.first*this->mRange*this->mRange + key.second*this->mRange + i;

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
