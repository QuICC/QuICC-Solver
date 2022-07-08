/** 
 * @file CurlNL.cpp
 * @brief Source of flag for curl of nonlinear term
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Path/CurlNL.hpp"

// Project includes
//
#include "QuICC/Transform/Path/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Path {

   std::string CurlNL::sTag()
   {
      return "Path::CurlNL";
   }

   const std::size_t& CurlNL::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<CurlNL>(CurlNL::sTag());
      return *i;
   }

   CurlNL::CurlNL()
      : IId(CurlNL::sTag())
   {
   }

   CurlNL::~CurlNL()
   {
   }

}
}
}
