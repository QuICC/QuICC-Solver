/**
 * @file Set3DModes.cpp
 * @brief Source ofspectral kernel to set 3D modes
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SpectralKernels/Set3DModes.hpp"

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Spectral {

namespace Kernel {

   Set3DModes::Set3DModes(const bool isComplex)
      : ISpectralKernel(isComplex)
   {
   }

   Set3DModes::~Set3DModes()
   {
   }

   void Set3DModes::init(const Real3DMapType modes)
   {
      if(this->mIsComplex)
      {
         throw std::logic_error("Initialized complex spectral kernel with real value.");
      }

      this->mspRModes = std::make_shared<Real3DMapType>();
      *this->mspRModes = modes;
   }

   void Set3DModes::init(const Complex3DMapType modes)
   {
      if(!this->mIsComplex)
      {
         throw std::logic_error("Initialized real spectral kernel with complex value.");
      }

      this->mspZModes = std::make_shared<Complex3DMapType>();
      *this->mspZModes = modes;
   }

   MHDVariant Set3DModes::compute(const int i, const int j, const int k) const
   {
      // Create key
      std::pair<int,int>   key;
      const auto& tRes = *this->spRes()->cpu()->dim(Dimensions::Transform::SPECTRAL);
      if(this->spRes()->sim().ss().has(SpatialScheme::Feature::SpectralOrdering132))
      {
         key = std::make_pair(tRes.idx<Dimensions::Data::DAT3D>(k),tRes.idx<Dimensions::Data::DAT2D>(j,k));
      } else
      {
         key = std::make_pair(tRes.idx<Dimensions::Data::DAT2D>(j,k),tRes.idx<Dimensions::Data::DAT3D>(k));
      }

      // Get value for key
      if(this->mIsComplex)
      {
         if(this->mspZModes->count(key) > 0)
         {
            if(this->mspZModes->find(key)->second.count(i) > 0)
            {
               return this->mspZModes->find(key)->second.find(i)->second;
            } else
            {
               return MHDComplex(0);
            }
         } else
         {
            return MHDComplex(0);
         }
      } else
      {
         if(this->mspRModes->count(key) > 0)
         {
            if(this->mspRModes->find(key)->second.count(i) > 0)
            {
               return this->mspRModes->find(key)->second.find(i)->second;
            } else
            {
               return MHDFloat(0);
            }
         } else
         {
            return MHDFloat(0);
         }
      }
   }

}
}
}
