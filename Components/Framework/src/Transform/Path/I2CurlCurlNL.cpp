/** 
 * @file I2CurlCurlNL.cpp
 * @brief Source of flag for second integral of curl curl of nonlinear term
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Path/I2CurlCurlNL.hpp"

// Project includes
//
#include "QuICC/Transform/Path/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Path {

   std::string I2CurlCurlNL::sTag()
   {
      return "Path::I2CurlCurlNL";
   }

   const std::size_t& I2CurlCurlNL::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I2CurlCurlNL>(I2CurlCurlNL::sTag());
      return *i;
   }

   I2CurlCurlNL::I2CurlCurlNL()
      : IId(I2CurlCurlNL::sTag())
   {
   }

   I2CurlCurlNL::~I2CurlCurlNL()
   {
   }

}
}
}
