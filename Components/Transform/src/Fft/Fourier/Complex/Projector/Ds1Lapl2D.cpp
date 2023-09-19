/**
 * @file Ds1Lapl2D.cpp
 * @brief Source of the implementation of the Fourier D(slow) 2D laplacian complex projector
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/Ds1Lapl2DBase.hpp"
#include "Types/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   void Ds1Lapl2D<base_t>::applyPreOperator(MatrixZ& tmp, const MatrixZ& in) const
   {
      std::vector<std::pair<int,int> > orders = { {2,1}, {0,3} };
      this->mBackend.inputDiff2D(tmp, in, orders, this->mspSetup->boxScale(), this->mspSetup->idBlocks());

   }

}
}
}
}
}
}
