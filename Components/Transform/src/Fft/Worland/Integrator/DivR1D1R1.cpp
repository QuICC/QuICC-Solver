/**
 * @file DivR1D1R1.cpp
 * @brief Source of the implementation of the Worland DivR1D1R1 integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Worland/Integrator/DivR1D1R1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   DivR1D1R1::DivR1D1R1()
   {
      this->setProfileTag();
   }

   DivR1D1R1::~DivR1D1R1()
   {
   }

   void DivR1D1R1::initBackend() const
   {
      int lshift = -1; // operator shifts l by one
      int extraN = 1; // 1 extra modes is required for l - 1 shift
      this->mBackend.init(*this->mspSetup, lshift, extraN);
      this->mBackend.addStorage(0, 1);
   }

   void DivR1D1R1::computeWorlandExpansion(const bool isEven) const
   {
      this->mBackend.forwardWorland(isEven);

      const int main = 0;
      const int extra = 1;

      // Copy expansion shifting indexes by -1 and scaling
      this->mBackend.copy(extra, main, -1, isEven);
      this->mBackend.scaleD(isEven, 0, extra);
      this->mBackend.lshift(extra, 1, isEven);
      this->mBackend.lowerAlpha(0.5, isEven, extra, 1.0);
      this->mBackend.lowerR2Beta(-0.5, isEven, extra, 1.0);
      this->mBackend.raiseR2Beta(-0.5, isEven, extra, 1.0/std::sqrt(2.0), true);

      // Scaling original expansion
      this->mBackend.scaleALPY(1.0, 1.0, isEven);
      this->mBackend.raiseR2Beta(-0.5, isEven, main);

      // Add second term
      this->mBackend.add(main, extra, 0, isEven);
   }

   void DivR1D1R1::applyPreOperator(const Matrix& in, const bool isEven) const
   {
      this->mBackend.input(in, isEven);
      this->mBackend.io(isEven);
   }

   void DivR1D1R1::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      this->computeWorlandExpansion(isEven);
      this->mBackend.output(rOut, isEven);
   }

   void DivR1D1R1::applyPreOperator(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      this->mBackend.input(in, isEven, useReal);
      this->mBackend.io(isEven);
   }

   void DivR1D1R1::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->computeWorlandExpansion(isEven);
      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
