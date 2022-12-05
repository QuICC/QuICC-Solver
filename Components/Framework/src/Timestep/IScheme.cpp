/**
 * @file IScheme.cpp
 * @brief Implementation of the interface to generic timestepping scheme
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Timestep/IScheme.hpp"

// Project includes
//

namespace QuICC {

namespace Timestep {

   IScheme::IScheme()
      : mUseEmbedded(false)
   {
   }

   // Use scheme's embedded lower order scheme?
   bool IScheme::useEmbedded() const
   {
      return this->mUseEmbedded;
   }

   void IScheme::enableEmbedded()
   {
      if(this->hasEmbedded())
      {
         this->mUseEmbedded = true;
      } else
      {
         throw std::logic_error("Tried to activate inexistant embedded scheme!");
      }
   }

}
}
