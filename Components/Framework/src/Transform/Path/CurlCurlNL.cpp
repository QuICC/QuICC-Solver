/** 
 * @file CurlCurlNL.cpp
 * @brief Source of flag for curl curl of nonlinear term
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Path/CurlCurlNL.hpp"

// Project includes
//
#include "QuICC/Transform/Path/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Path {

   std::string CurlCurlNL::sTag()
   {
      return "Path::CurlCurlNL";
   }

   const std::size_t& CurlCurlNL::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<CurlCurlNL>(CurlCurlNL::sTag());
      return *i;
   }

   CurlCurlNL::CurlCurlNL()
      : IId(CurlCurlNL::sTag())
   {
   }

   CurlCurlNL::~CurlCurlNL()
   {
   }

}
}
}
