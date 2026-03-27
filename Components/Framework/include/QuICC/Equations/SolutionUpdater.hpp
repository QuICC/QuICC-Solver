/**
 * @file SolutionUpdater.hpp
 * @brief Basic solution updater (simple passthrough)
 */

#pragma once

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Equations {

/**
 * @brief Basic solution updater (passthrough)
 */
class SolutionUpdater
{
   public:
      /**
       * @brief ctor
       */
      SolutionUpdater() = default;

      /**
       * @brief dtor
       */
      virtual ~SolutionUpdater() = default;

      /**
       * @brief Simple passthrough
       */
      MHDVariant operator()(const MHDVariant newData, const int, const int, const int) const {return newData;};
};

} // namespace Equations
} // namespace QuICC
