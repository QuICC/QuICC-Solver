/**
 * @file R1.cpp
 * @brief Source of the implementation of the Worland R1 integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Worland/Integrator/R1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   R1::R1()
   {
      this->mProfileId = Debug::Profiler::WORLANDINTG_R1;
   }

   R1::~R1()
   {
   }

   void R1::initBackend() const
   {
      int lshift = -1; // operator shifts l by -1
      int extraN = 1; // 1 extra modes is required
      this->mBackend.init(*this->mspSetup, lshift, extraN);
   }

   void R1::computeWorlandExpansion(const bool isEven) const
   {
      this->mBackend.forwardWorland(isEven);
      this->mBackend.raiseBeta(-0.5, isEven);
   }

   void R1::applyPreOperator(const Matrix& in, const bool isEven) const
   {
      this->mBackend.input(in, isEven);
      this->mBackend.io(isEven);
   }

   void R1::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      this->computeWorlandExpansion(isEven);
      this->mBackend.output(rOut, isEven);
   }

   void R1::applyPreOperator(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      this->mBackend.input(in, isEven, useReal);
      this->mBackend.io(isEven);
   }

   void R1::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->computeWorlandExpansion(isEven);
      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
