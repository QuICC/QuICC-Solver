/** 
 * @file Sin_1LaplhDphi.cpp
 * @brief Source of backward projection operator Sin_1LaplhDphi 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/Sin_1LaplhDphi.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string Sin_1LaplhDphi::sTag()
   {
      return "Bwd::Sin_1LaplhDphi";
   }

   const std::size_t& Sin_1LaplhDphi::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Sin_1LaplhDphi>(Sin_1LaplhDphi::sTag());
      return *i;
   }

   Sin_1LaplhDphi::Sin_1LaplhDphi()
      : IOperator(Sin_1LaplhDphi::sTag())
   {
   }

   Sin_1LaplhDphi::~Sin_1LaplhDphi()
   {
   }

}
}
}
