/** 
 * @file I4CurlCurlNL.cpp
 * @brief Source of flag for fourth integral of curl curl of nonlinear term
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Path/I4CurlCurlNL.hpp"

// Project includes
//
#include "QuICC/Transform/Path/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Path {

   std::string I4CurlCurlNL::sTag()
   {
      return "Path::I4CurlCurlNL";
   }

   const std::size_t& I4CurlCurlNL::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I4CurlCurlNL>(I4CurlCurlNL::sTag());
      return *i;
   }

   I4CurlCurlNL::I4CurlCurlNL()
      : IId(I4CurlCurlNL::sTag())
   {
   }

   I4CurlCurlNL::~I4CurlCurlNL()
   {
   }

}
}
}
