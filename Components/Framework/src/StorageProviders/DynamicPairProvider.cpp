/** 
 * @file DynamicPairProvider.cpp
 * @brief Source of a data pair storage provider adapting its size dynamically
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/StorageProviders/DynamicPairProvider.hpp"

// Project includes
//

namespace QuICC {

   DynamicPairProvider::DynamicPairProvider()
   {
   }

   DynamicPairProvider::~DynamicPairProvider()
   {
   }

   void DynamicPairProvider::init(SharedFwdSetupType spSetupFwd, SharedBwdSetupType spSetupBwd)
   {
      this->mFwd.init(spSetupFwd);

      this->mBwd.init(spSetupBwd);
   }

   MHDFloat  DynamicPairProvider::requiredStorage() const
   {
      MHDFloat mem = 0.0;

   #ifdef QUICC_STORAGEPROFILE
      mem += this->mFwd.requiredStorage();

      mem += this->mBwd.requiredStorage();
   #endif // QUICC_STORAGEPROFILE

      return mem;
   }


}
