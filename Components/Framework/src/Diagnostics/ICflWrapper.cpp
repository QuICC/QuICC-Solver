/**
 * @file ICflWrapper.cpp
 * @brief Source of the interface for the CFL constraint
 */

// Debug includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Diagnostics/ICflWrapper.hpp"

// Project includes
//

namespace QuICC {

namespace Diagnostics {

   ICflWrapper::ICflWrapper(const SharedIVectorWrapper spVelocity)
      : mspVelocity(spVelocity), mspMagnetic()
   {
   }

   ICflWrapper::ICflWrapper(const SharedIVectorWrapper spVelocity, const SharedIVectorWrapper spMagnetic)
      : mspVelocity(spVelocity), mspMagnetic(spMagnetic)
   {
   }

   ICflWrapper::~ICflWrapper()
   {
   }

   void ICflWrapper::updateCflMatrix(Matrix& cfl) const
   {
      int idx;
      cfl(0,0) = cfl.row(0).tail(cfl.cols()-1).minCoeff(&idx);

      if(cfl.rows() > 1)
      {
         cfl.col(0).tail(cfl.rows()-1) = cfl.col(idx+1).tail(cfl.rows()-1);
      }
   }

}
}
