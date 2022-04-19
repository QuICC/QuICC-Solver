/**
 * @file P.cpp
 * @brief Source of the implementation of the Worland P projector
 */

// System includes
//
#include <cassert>

// Debug includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Worland/Projector/P.hpp"

// Project includes
//
#include "QuICC/Debug/Profiler/ProfilerMacro.h"
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Projector {

   P::P()
   {
      this->mProfileId = Debug::Profiler::WORLANDPROJ_P;
   }

   P::~P()
   {
   }

   void P::computeWorlandExpansion(const bool isEven) const
   {
      this->mBackend.scaleC(1.0/std::sqrt(Math::PI), isEven);
      this->mBackend.backwardWorland(isEven);
   }

   void P::applyPreOperator(const Matrix& in, const bool isEven) const
   {
      this->mBackend.input(in, isEven, true);
      this->computeWorlandExpansion(isEven);
      this->mBackend.io(isEven);
   }

   void P::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      this->mBackend.output(rOut, isEven);
   }

   void P::applyPreOperator(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      this->mBackend.input(in, isEven, useReal, true);
      this->computeWorlandExpansion(isEven);
      this->mBackend.io(isEven);
   }

   void P::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
