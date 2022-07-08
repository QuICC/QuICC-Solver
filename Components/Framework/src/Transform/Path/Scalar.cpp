/** 
 * @file Scalar.cpp
 * @brief Source of flag for scalar of physical values
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Path/Scalar.hpp"

// Project includes
//
#include "QuICC/Transform/Path/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Path {

   std::string Scalar::sTag()
   {
      return "Path::Scalar";
   }

   const std::size_t& Scalar::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Scalar>(Scalar::sTag());
      return *i;
   }

   Scalar::Scalar()
      : IId(Scalar::sTag())
   {
   }

   Scalar::~Scalar()
   {
   }

}
}
}
