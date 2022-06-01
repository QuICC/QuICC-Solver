/** 
 * @file IEnvironment.cpp
 * @brief Source of generic environment
 */

// Configuration includes
//
#include "QuICC/Debug/DebuggerMacro.h"

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Environment/IEnvironment.hpp"

// Project includes
//

namespace QuICC {

namespace Environment {

   int IEnvironment::mIoRank = -99;

   int IEnvironment::mSize = -99;

   int IEnvironment::mId = -99;

   IEnvironment::IEnvironment()
   {
   }

   IEnvironment::~IEnvironment()
   {
   }

   int IEnvironment::ioRank() const
   {
      return this->mIoRank;
   }

   bool IEnvironment::allowsIO() const
   {
      return (this->ioRank() == this->id());
   }

   void IEnvironment::checkEnvironment(const int size)
   {
      // Check that the size was set
      if(this->size() < 1)
      {
         this->abort("Environment contains has negative size!");
      }

      // Check that the local ID was set
      if(this->id() < 0)
      {
         this->abort("Environment local ID was not set!");
      }

      // Check that IO rank was set
      if(this->ioRank() < 0)
      {
         this->abort("Environment IO rank was not set!");
      }

      // Check compatibility between requested cores and setup cores
      if(size != this->size())
      {
         this->abort("Environment parameters and setup have conflicting sizes: "
            +std::to_string(size)+" vs "+std::to_string(this->size()));
      }
   }

}
}
