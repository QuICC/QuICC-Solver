/** 
 * @file DynamicRZProvider.cpp
 * @brief Source of a combined real and complex data storage provider adapting its size dynamically
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/StorageProviders/DynamicRZProvider.hpp"

// Project includes
//

namespace QuICC {

   DynamicRZProvider::DynamicRZProvider()
   {
   }

   DynamicRZProvider::~DynamicRZProvider()
   {
   }

   void DynamicRZProvider::init(SharedSetupType spSetup)
   {
      // initialize real storage
      this->mReal.init(spSetup);

      // initialize complex storage
      this->mComplex.init(spSetup);
   }

   MHDFloat  DynamicRZProvider::requiredStorage() const
   {
      MHDFloat mem = 0.0;
      
   #ifdef QUICC_STORAGEPROFILE
      mem += this->mReal.requiredStorage();

      mem += this->mComplex.requiredStorage();
   #endif // QUICC_STORAGEPROFILE

      return mem;
   }

}
