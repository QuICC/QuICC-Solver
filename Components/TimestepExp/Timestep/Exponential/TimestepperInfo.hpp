/**
 * @file SparseCoordinatorData.hpp
 * @brief Implementation of the base for a general sparse solver coordinator
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_TIMESTEPPERINFO_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_TIMESTEPPERINFO_HPP

// System includes
//
#include <memory>
#include <map>

// Project includes
//
#include "Types/Typedefs.hpp"

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
      std::size_t matStart;
      std::vector<std::size_t> matIds;
      std::map<std::size_t, std::map<std::size_t, std::pair<int, DecoupledZSparse>>> ops;
   };

} // Exponential
} // Timestep
} // QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_TIMESTEPPERINFO_HPP
