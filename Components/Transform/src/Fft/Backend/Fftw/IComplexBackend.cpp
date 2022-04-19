/**
 * @file IComplexBackend.cpp
 * @brief Source of the interface for a generic FFTW based complex projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/Fftw/IComplexBackend.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   IComplexBackend::IComplexBackend()
   {
   }

   IComplexBackend::~IComplexBackend()
   {
   }

   void IComplexBackend::init(const SetupType& setup) const
   {
      // Get size of positive and negative frequency parts
      int specSize = setup.specSize();
      this->mNegN = specSize/2;
      this->mPosN = this->mNegN + (specSize%2);
   }

   void IComplexBackend::initMeanBlocks(const MatrixI& idBlocks) const
   {
      // Set the mean blocks
      int start = 0;
      for(int i = 0; i < idBlocks.rows(); ++i)
      {
         if(idBlocks(i,0) == 0)
         {
            this->mMeanBlocks.push_back(std::make_pair(start, idBlocks(i,1)));
         }

         start += idBlocks(i,1);
      }
   }

   void IComplexBackend::input(const MHDComplex* in) const
   {
      this->io(this->mTmp.data(), in);
   }

   void IComplexBackend::output(MHDComplex* out) const
   {
      this->io(out, this->mTmp.data());
   }

   void IComplexBackend::io(MHDComplex* out, const MHDComplex* in) const
   {
      this->mpOut = out;
      this->mpIn = in;
   }

   Array IComplexBackend::positiveK() const
   {
      return Array::LinSpaced(this->mPosN, 0, this->mPosN-1);
   }

   Array IComplexBackend::negativeK() const
   {
      return Array::LinSpaced(this->mNegN, 0, this->mNegN-1).array() - static_cast<MHDFloat>(this->mNegN);
   }

}
}
}
}
}
