/**
 * @file Lapl2D.cpp
 * @brief Source of the implementation of the Fourier complex D projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/Lapl2D.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   Lapl2D::Lapl2D()
   {
   }

   Lapl2D::~Lapl2D()
   {
   }

   void Lapl2D::applyPreOperator(MatrixZ& rOut, const MatrixZ& in) const
   {
      std::vector<std::pair<int,int> > orders = { {2,0}, {0,2} };
      this->mBackend.inputDiff2D(in, orders, this->mspSetup->boxScale(), this->mspSetup->idBlocks());

      this->mBackend.output(rOut.data());
   }

   void Lapl2D::applyPostOperator(MatrixZ&) const
   {
   }

}
}
}
}
}
}
