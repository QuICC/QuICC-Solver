/**
 * @file VariableBase.cpp
 * @brief Base of the implementation of the variables
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Variables/VariableBase.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Datatypes {

   VariableBase::VariableBase(SharedResolution spRes)
      : mspRes(spRes)
   {
   }

   VariableBase::~VariableBase()
   {
   }

   MHDFloat VariableBase::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

}
}
