/**
 * @file MinEqInfo.hpp
 * @brief 
 */

#pragma once

// System includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

/**
 * @brief Minimal equation information
 */
struct MinEqInfo
{
   /**
    * @brief ctor
    */
   MinEqInfo(const int it, const std::size_t timeId, const SpectralFieldId& myId): it(it), timeId(timeId), myId(myId){};

   /**
    * @brief dtor
    */
   ~MinEqInfo() = default;

   /// Iteration
   int it;

   /// Active solver time
   std::size_t timeId;

   /// Field ID
   SpectralFieldId myId;
};

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
