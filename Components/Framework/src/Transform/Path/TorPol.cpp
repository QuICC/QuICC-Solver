/** 
 * @file TorPol.cpp
 * @brief Source of flag for Toroidal/Poloidal of physical values
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Path/TorPol.hpp"

// Project includes
//
#include "QuICC/Transform/Path/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Path {

   std::string TorPol::sTag()
   {
      return "Path::TorPol";
   }

   const std::size_t& TorPol::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<TorPol>(TorPol::sTag());
      return *i;
   }

   TorPol::TorPol()
      : IId(TorPol::sTag())
   {
   }

   TorPol::~TorPol()
   {
   }

}
}
}
