/** 
 * @file I2CurlNL.cpp
 * @brief Source of flag for second integral of curl of nonlinear term
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Path/I2CurlNL.hpp"

// Project includes
//
#include "QuICC/Transform/Path/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Path {

   std::string I2CurlNL::sTag()
   {
      return "Path::I2CurlNL";
   }

   const std::size_t& I2CurlNL::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I2CurlNL>(I2CurlNL::sTag());
      return *i;
   }

   I2CurlNL::I2CurlNL()
      : IId(I2CurlNL::sTag())
   {
   }

   I2CurlNL::~I2CurlNL()
   {
   }

}
}
}
