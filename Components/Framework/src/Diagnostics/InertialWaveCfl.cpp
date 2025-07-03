/**
 * @file InertialWaveCfl.cpp
 * @brief Source of the CFL constraint wrapper in a spherical geometry
 */

// System includes
//

// Project includes
//
#include "QuICC/Diagnostics/InertialWaveCfl.hpp"
#include "QuICC/NonDimensional/CflInertial.hpp"

namespace QuICC {

namespace Diagnostics {

   InertialWaveCfl::InertialWaveCfl(const std::map<std::size_t,NonDimensional::SharedINumber>& params, const MHDFloat courant)
      : ICflWrapper(courant), mCfl(2,1)
   {
      // Inertial wave CFL
      if(params.count(NonDimensional::CflInertial::id()) > 0)
      {
         this->mCfl(0,0) = params.find(NonDimensional::CflInertial::id())->second->value();
         this->mIsActive = true;
      }
      else
      {
         this->mCfl(0,0) = -1;
         this->mIsActive = false;
      }
      this->mCfl(1,0) = -1;
   }

   void InertialWaveCfl::init(const std::vector<Array>& mesh)
   {
   }

   Matrix InertialWaveCfl::initialCfl() const
   {
      return this->mCfl;
   }

   Matrix InertialWaveCfl::cfl() const
   {
      return this->mCfl;
   }

}
}
