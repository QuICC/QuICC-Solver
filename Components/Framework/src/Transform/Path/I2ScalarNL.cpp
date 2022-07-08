/** 
 * @file I2ScalarNL.cpp
 * @brief Source of flag for scalar of nonlinear term
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Path/I2ScalarNL.hpp"

// Project includes
//
#include "QuICC/Transform/Path/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Path {

   std::string I2ScalarNL::sTag()
   {
      return "Path::I2ScalarNL";
   }

   const std::size_t& I2ScalarNL::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I2ScalarNL>(I2ScalarNL::sTag());
      return *i;
   }

   I2ScalarNL::I2ScalarNL()
      : IId(I2ScalarNL::sTag())
   {
   }

   I2ScalarNL::~I2ScalarNL()
   {
   }

}
}
}
