/**
 * @file Add.hpp
 * @brief Add arithmetics
 */

#ifndef QUICC_ARITHMETICS_ADD_HPP
#define QUICC_ARITHMETICS_ADD_HPP

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
    * @brief Add Arithmetics
    */
   class Add: public IRegisterId<Add>
   {
      public:
         /**
          * @brief Constructor
          */
         Add();

         friend class IRegisterId<Add>;

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

#endif // QUICC_ARITHMETICS_ADD_HPP
