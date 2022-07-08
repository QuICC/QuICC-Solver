/**
 * @file DivR1.cpp
 * @brief Source of the implementation of the Worland DivR1 integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Worland/Integrator/DivR1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   DivR1::DivR1()
   {
      this->mProfileId = Debug::Profiler::WORLANDINTG_DIVR1;
   }

   DivR1::~DivR1()
   {
   }

   void DivR1::initBackend() const
   {
      int lshift = -1; // operator shifts l by one
      int extraN = 0; // no extra modes are required
      this->mBackend.init(*this->mspSetup, lshift, extraN);
   }

   void DivR1::computeWorlandExpansion(const bool isEven) const
   {
      this->mBackend.forwardWorland(isEven);
      this->mBackend.raiseR2Beta(-0.5, isEven);
   }

   void DivR1::applyPreOperator(const Matrix& in, const bool isEven) const
   {
      this->mBackend.input(in, isEven);
      this->mBackend.io(isEven);
   }

   void DivR1::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      this->computeWorlandExpansion(isEven);
      this->mBackend.output(rOut, isEven);
   }

   void DivR1::applyPreOperator(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      this->mBackend.input(in, isEven, useReal);
      this->mBackend.io(isEven);
   }

   void DivR1::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->computeWorlandExpansion(isEven);
      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
