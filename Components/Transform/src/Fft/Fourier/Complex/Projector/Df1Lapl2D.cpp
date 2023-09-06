/**
 * @file Df1Lapl2D.cpp
 * @brief Source of the implementation of the Fourier  D(fast) 2D laplacian complex projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/Df1Lapl2D.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   void Df1Lapl2D<base_t>::applyPreOperator(MatrixZ& tmp, const MatrixZ& in) const
   {
      std::vector<std::pair<int,int> > orders = { {3,0}, {1,2} };
      this->mBackend.inputDiff2D(tmp, in, orders, this->mspSetup->boxScale(), this->mspSetup->idBlocks());
   }

}
}
}
}
}
}
