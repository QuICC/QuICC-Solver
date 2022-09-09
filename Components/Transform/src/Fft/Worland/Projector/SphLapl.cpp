/**
 * @file SphLapl.cpp
 * @brief Source of the implementation of the Worland SphLapl projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Worland/Projector/SphLapl.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Projector {

   SphLapl::SphLapl()
   {
      this->setProfileTag();
   }

   SphLapl::~SphLapl()
   {
   }

   void SphLapl::initBackend() const
   {
      int lshift = 0; // operator doesn't shift l
      int extraN = 0; // no extra modes are required
      this->mBackend.init(*this->mspSetup, lshift, extraN);
      this->mBackend.addStorage(1, 0);
   }

   void SphLapl::computeWorlandExpansion(const bool isEven) const
   {
      const int main = 0;
      const int extra = 1;

      this->mBackend.scaleC(1.0/std::sqrt(Math::PI), isEven);

      // Copy expansion shifting indexes by -2 and scaling
      this->mBackend.copy(extra, main, -2, isEven);
      this->mBackend.scaleSphLaplA(isEven, 0, extra);
      this->mBackend.lshift(extra, 2, isEven);
      this->mBackend.lowerAlpha(1.5, isEven, extra, 1.0);
      this->mBackend.lowerR2Beta(0.5, isEven, extra, 1.0);

      // Shifting original expansion by -1 and scaling
      this->mBackend.nshift(main, -1, isEven);
      this->mBackend.scaleSphLaplB(isEven);

      // Add second term
      this->mBackend.lshift(main, 1, isEven);
      this->mBackend.add(main, extra, 0, isEven);
      this->mBackend.lowerAlpha(0.5, isEven, main);
      this->mBackend.lowerBeta(-0.5, isEven, main);

      this->mBackend.backwardWorland(isEven);
   }

   void SphLapl::applyPreOperator(const Matrix& in, const bool isEven) const
   {
      this->mBackend.input(in, isEven, true);
      this->computeWorlandExpansion(isEven);
      this->mBackend.io(isEven);
   }

   void SphLapl::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      this->mBackend.output(rOut, isEven);
   }

   void SphLapl::applyPreOperator(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      this->mBackend.input(in, isEven, useReal, true);
      this->computeWorlandExpansion(isEven);
      this->mBackend.io(isEven);
   }

   void SphLapl::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
