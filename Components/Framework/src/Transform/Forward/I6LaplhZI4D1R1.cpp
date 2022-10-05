/**
 * @file I6LaplhZI4D1R1.cpp
 * @brief Source of the forward transform operator Forward::I6LaplhZI4D1R1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I6LaplhZI4D1R1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I6LaplhZI4D1R1::sTag()
   {
      return "Fwd::I6LaplhZI4D1R1";
   }

   std::string I6LaplhZI4D1R1::sFormatted()
   {
      return "Forward::I6LaplhZI4D1R1";
   }

   I6LaplhZI4D1R1::I6LaplhZI4D1R1()
      : IRegisterId<I6LaplhZI4D1R1>(I6LaplhZI4D1R1::sTag(), I6LaplhZI4D1R1::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
