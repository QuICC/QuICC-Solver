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
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/Lapl2DBase.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   void Lapl2D<base_t>::applyPreOperator(MatrixZ& tmp, const MatrixZ& in) const
   {
      std::vector<std::pair<int,int> > orders = { {2,0}, {0,2} };
      this->mBackend.inputDiff2D(tmp, in, orders, this->mspSetup->boxScale(), this->mspSetup->idBlocks());

   }

}
}
}
}
}
}
