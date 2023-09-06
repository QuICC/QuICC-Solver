/**
 * @file InvLapl2D.cpp
 * @brief Source of the implementation of the Fourier complex inverse 2D Laplacian integrator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/InvLapl2DBase.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {

   void InvLapl2D<base_t>::applyPostOperator(MatrixZ& rOut) const
   {
      std::vector<std::pair<int,int> > orders = { {2,0}, {0,2} };
      int id = this->mBackend.computeDiff2D(orders, this->mspSetup->boxScale(), this->mspSetup->idBlocks(), true);
      this->mBackend.applyDiff2D(rOut, id);
      this->mBackend.destroyDiff2D(id);
   }

} // namespace Integrator
} // namespace Complex
} // namespace Fourier
} // namespace Fft
} // namespace Transform
} // namespace QuICC
