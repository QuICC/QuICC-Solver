/**
 * @file Wnl.cpp
 * @brief Source of the implementation of Wnl
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Polynomial/Worland/Wnl.hpp"

// Project includes
//

namespace QuICC {

namespace Polynomial {

namespace Worland {

   Wnl::Wnl(const Internal::MHDFloat alpha, const Internal::MHDFloat dBeta, const int lShift)
      : WorlandBase(alpha, dBeta), mLShift(lShift)
   {
   }

}
}
}

