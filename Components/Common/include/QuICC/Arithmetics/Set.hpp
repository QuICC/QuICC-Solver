/**
 * @file Set.hpp
 * @brief Set arithmetics
 */

#ifndef QUICC_ARITHMETICS_SET_HPP
#define QUICC_ARITHMETICS_SET_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Arithmetics/IRegisterId.hpp"

namespace QuICC {

namespace Arithmetics {

   /**
    * @brief Set Arithmetics
    */
   class Set: public IRegisterId<Set>
   {
      public:
         /**
          * @brief Constructor
          */
         Set();

         friend class IRegisterId<Set>;

      protected:

      private:
         /**
          * @brief Unique tag
          */
         static std::string sTag();

         /**
          * @brief Formatted name
          */
         static std::string sFormatted();
   };

}
}

#endif // QUICC_ARITHMETICS_SET_HPP
