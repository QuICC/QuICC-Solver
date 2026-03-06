/**
 * @file SparseCoordinatorData.hpp
 * @brief Implementation of the base for a general sparse solver coordinator
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_TIMESTEPPERINFO_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_TIMESTEPPERINFO_HPP

// System includes
//
#include <memory>

// Project includes
//

namespace QuICC {

namespace Timestep {

namespace Exponential {

   struct TimestepperInfo
   {
      bool isComplex;
      std::size_t solverIndex;
      std::size_t fieldIndex;
      std::size_t rows;
      std::size_t cols;
      std::size_t blockN;
   };

} // Exponential
} // Timestep
} // QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_TIMESTEPPERINFO_HPP
