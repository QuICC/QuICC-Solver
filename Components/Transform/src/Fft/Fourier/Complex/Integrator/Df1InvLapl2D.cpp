/**
 * @file Df1InvLapl2D.cpp
 * @brief Source of the implementation of the Fourier complex D(fast) of inverse 2D Laplacian integrator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/Df1InvLapl2DBase.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {

   void Df1InvLapl2D<base_t>::applyPostOperator(MatrixZ& rOut) const
   {
      std::vector<std::pair<int,int> > orders = { {2,0}, {0,2} };
      int invId = this->mBackend.computeDiff2D(orders, this->mspSetup->boxScale(), this->mspSetup->idBlocks(), true);
      orders.clear();
      orders = { {1,0} };
      int dfId = this->mBackend.computeDiff2D(orders, this->mspSetup->boxScale(), this->mspSetup->idBlocks());
      int opId = this->mBackend.multDiff2D(invId, dfId);
      this->mBackend.destroyDiff2D(invId);
      this->mBackend.destroyDiff2D(dfId);
      this->mBackend.applyDiff2D(rOut, opId);
      this->mBackend.destroyDiff2D(opId);
   }

} // namespace Integrator
} // namespace Complex
} // namespace Fourier
} // namespace Fft
} // namespace Transform
} // namespace QuICC
