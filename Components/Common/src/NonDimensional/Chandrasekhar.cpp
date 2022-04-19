/**
 * @file Chandrasekhar.cpp
 * @brief Source of the Chandrasekhar nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Chandrasekhar.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Chandrasekhar::sTag()
   {
      return "chandrasekhar";
   }

   std::string Chandrasekhar::sFormatted()
   {
      return "Chandrasekhar";
   }

   Chandrasekhar::Chandrasekhar(const MHDFloat value)
      : IRegisterId<Chandrasekhar>(value, Chandrasekhar::sTag(), Chandrasekhar::sFormatted())
   {
   }

}
}
