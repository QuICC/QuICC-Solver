/**
 * @file ICflWrapper.cpp
 * @brief Source of the interface for the CFL constraint
 */

// System includes
//

// Project includes
//
#include "QuICC/Diagnostics/ICflWrapper.hpp"

namespace QuICC {

namespace Diagnostics {

   ICflWrapper::ICflWrapper(const MHDFloat courant)
      : mcCourant(courant), mIsActive(false)
   {
   }

   bool ICflWrapper::isActive() const
   {
      return this->mIsActive;
   }

   void ICflWrapper::setField(const std::size_t id, const SharedIVectorWrapper spField)
   {
      if(this->mFields.count(id) != 1)
      {
         throw std::logic_error("Unexpected field given to CFL wrapper");
      }

      this->mFields.at(id) = spField;
   }

} // namespace Diagnostics
} // namespace QuICC
