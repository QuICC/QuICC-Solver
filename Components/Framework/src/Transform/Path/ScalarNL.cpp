/** 
 * @file ScalarNL.cpp
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
#include "QuICC/Transform/Path/ScalarNL.hpp"

// Project includes
//
#include "QuICC/Transform/Path/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Path {

   std::string ScalarNL::sTag()
   {
      return "Path::ScalarNL";
   }

   const std::size_t& ScalarNL::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<ScalarNL>(ScalarNL::sTag());
      return *i;
   }

   ScalarNL::ScalarNL()
      : IId(ScalarNL::sTag())
   {
   }

   ScalarNL::~ScalarNL()
   {
   }

}
}
}
