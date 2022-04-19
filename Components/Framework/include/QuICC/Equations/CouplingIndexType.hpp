/**
 * @file CouplingIndexType.hpp
 * @brief Enum to specify the type of the indexes
 */

#ifndef QUICC_EQUATIONS_COUPLINGINDEXTYPE_HPP
#define QUICC_EQUATIONS_COUPLINGINDEXTYPE_HPP

// Configuration includes
//

// System includes
//
#include <set>
#include <stdexcept>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Equations {

   /**
    * @brief Enum to specify the type of indexes
    */
   enum class CouplingIndexType {
      /// Matrix index is slowest index of field
      SLOWEST_SINGLE_RHS = 0,
      /// Matrix index is slowest index of field but has multiple RHS
      SLOWEST_MULTI_RHS,
      /// Matrix index is a mode index
      MODE,
      /// Single matrix (ex. TTT scheme)
      SINGLE,
   };

   inline CouplingIndexType safe_CouplingIndexType_cast(const std::size_t id);

   CouplingIndexType safe_CouplingIndexType_cast(const std::size_t id)
   {
      if(id == static_cast<int>(CouplingIndexType::SLOWEST_SINGLE_RHS) ||
         id == static_cast<int>(CouplingIndexType::SLOWEST_MULTI_RHS) ||
         id == static_cast<int>(CouplingIndexType::MODE) ||
         id == static_cast<int>(CouplingIndexType::SINGLE))
      {
         return static_cast<CouplingIndexType>(id);
      } else
      {
         throw std::logic_error("integer out of range for CouplingIndexType");
      }
   }

}
}

#endif // QUICC_EQUATIONS_COUPLINGINDEXTYPE_HPP
