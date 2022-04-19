/**
 * @file Library.cu
 * @brief Source for the static interface to cuFFT
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/Library.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   int Library::sCounter = 0;

   Library::Library()
   {
   }

   Library::~Library()
   {
   }

   void Library::registerUser()
   {
      ++Library::sCounter;
   }

   void Library::unregisterUser()
   {
      --Library::sCounter;
   }

   void Library::init()
   {
   }

   void Library::cleanup()
   {
      // Check if all objects have been destroyed
      if(Library::sCounter == 0)
      {
      }
   }

}
}
}
}
}
