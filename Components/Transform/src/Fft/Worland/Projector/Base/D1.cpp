/**
 * @file D1.cpp
 * @brief Source of the implementation of the Worland D1 projector
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Projector/Base/D1.hpp"
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Projector {

   D1<base_t>::D1()
   {
      this->setProfileTag();
   }

   void D1<base_t>::initBackend() const
   {
      int lshift = -1; // operator shifts l by -1
      int extraN = 1; // 1 extra mode is required
      bool onlyShiftParity = true;
      this->mBackend.init(*this->mspSetup, lshift, extraN, onlyShiftParity);
      this->mBackend.addStorage(1, 0);
   }

   void D1<base_t>::computeWorlandExpansion(const bool isEven) const
   {
      const int main = 0;
      const int extra = 1;

      this->mBackend.scaleC(1.0/std::sqrt(Math::PI), isEven);

      // Copy expansion shifting indexes by -1 and scaling
      this->mBackend.copy(extra, main, -1, isEven);
      this->mBackend.scaleD(isEven, 0, extra);
      this->mBackend.lshift(extra, 1, isEven);
      this->mBackend.lowerAlpha(0.5, isEven, extra, 1.0);
      this->mBackend.lowerR2Beta(-0.5, isEven, extra, 1.0);

      // Scaling original expansion
      this->mBackend.scaleALPY(1.0, 0.0, isEven);

      // Add second term
      this->mBackend.add(main, extra, 0, isEven);
      this->mBackend.lowerBeta(-0.5, isEven);

      this->mBackend.backwardWorland(isEven);
   }

   void D1<base_t>::applyPreOperator(const Matrix& in, const bool isEven) const
   {
      this->mBackend.input(in, isEven, true);
      this->computeWorlandExpansion(isEven);
      this->mBackend.io(isEven);
   }

   void D1<base_t>::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      this->mBackend.output(rOut, isEven);
   }

   void D1<base_t>::applyPreOperator(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      this->mBackend.input(in, isEven, useReal, true);
      this->computeWorlandExpansion(isEven);
      this->mBackend.io(isEven);
   }

   void D1<base_t>::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
