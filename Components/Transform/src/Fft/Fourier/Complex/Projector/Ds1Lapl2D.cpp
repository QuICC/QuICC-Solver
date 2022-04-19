/**
 * @file Ds1Lapl2D.cpp
 * @brief Source of the implementation of the Fourier D(slow) 2D laplacian complex projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/Ds1Lapl2D.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   Ds1Lapl2D::Ds1Lapl2D()
   {
   }

   Ds1Lapl2D::~Ds1Lapl2D()
   {
   }

   void Ds1Lapl2D::applyPreOperator(MatrixZ& rOut, const MatrixZ& in) const
   {
      std::vector<std::pair<int,int> > orders = { {2,1}, {0,3} };
      this->mBackend.inputDiff2D(in, orders, this->mspSetup->boxScale(), this->mspSetup->idBlocks());

      this->mBackend.output(rOut.data());
   }

   void Ds1Lapl2D::applyPostOperator(MatrixZ&) const
   {
   }

}
}
}
}
}
}
