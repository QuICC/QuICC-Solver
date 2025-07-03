/**
 * @file TorsionalOscillationCfl.cpp
 * @brief Source of the CFL constraint wrapper in a spherical geometry
 */

// System includes
//

// Project includes
//
#include "QuICC/Diagnostics/TorsionalOscillationCfl.hpp"
#include "QuICC/NonDimensional/CflTorsional.hpp"

namespace QuICC {

namespace Diagnostics {

   TorsionalOscillationCfl::TorsionalOscillationCfl(const std::map<std::size_t,NonDimensional::SharedINumber>& params, const MHDFloat courant)
      : ICflWrapper(courant), mCfl(2,1)
   {
      // Torsional wave CFL
      if(params.count(NonDimensional::CflTorsional::id()) > 0)
      {
         this->mCfl(0,0) = params.find(NonDimensional::CflTorsional::id())->second->value();
         this->mIsActive = false;
      }
      else
      {
         this->mCfl(0,0) = -1;
         this->mIsActive = false;
      }
      this->mCfl(1,0) = -1;

   }

   void TorsionalOscillationCfl::init(const std::vector<Array>& mesh)
   {
   }

   Matrix TorsionalOscillationCfl::initialCfl() const
   {
      return this->mCfl;
   }

   Matrix TorsionalOscillationCfl::cfl() const
   {
      return this->mCfl;
   }

}
}
