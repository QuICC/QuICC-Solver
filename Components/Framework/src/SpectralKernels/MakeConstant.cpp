/**
 * @file MakeConstant.cpp
 * @brief Source of trivial spectral kernel to make field constant
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SpectralKernels/MakeConstant.hpp"

// Project includes
//

namespace QuICC {

namespace Spectral {

namespace Kernel {

   MakeConstant::MakeConstant(const bool isComplex)
      : ISpectralKernel(isComplex)
   {
   }

   MakeConstant::~MakeConstant()
   {
   }

   void MakeConstant::init(const MHDFloat value)
   {
      if(this->mIsComplex)
      {
         throw std::logic_error("Initialized complex spectral kernel with real value.");
      }

      this->mRValue = value;
   }

   void MakeConstant::init(const MHDComplex value)
   {
      if(!this->mIsComplex)
      {
         throw std::logic_error("Initialized real spectral kernel with complex value.");
      }

      this->mZValue = value;
   }

   MHDVariant MakeConstant::compute(const int i, const int j, const int k) const
   {
      if(this->mIsComplex)
      {
         return this->mZValue;
      } else
      {
         return this->mRValue;
      }
   }

}
}
}
