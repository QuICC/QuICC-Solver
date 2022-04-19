/**
 * @file ITransform.cpp
 * @brief Source of the generic transform interface
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/ITransform.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

   void ITransform::unimplemented()
   {
      throw std::logic_error("Transform does not implemnent this option");
   }
}
}
