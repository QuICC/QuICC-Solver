/**
 * @file SetNeg.hpp
 * @brief SetNeg arithmetics
 */

#ifndef QUICC_ARITHMETICS_SETNEG_HPP
#define QUICC_ARITHMETICS_SETNEG_HPP

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
    * @brief SetNeg Arithmetics
    */
   class SetNeg: public IRegisterId<SetNeg>
   {
      public:
         /**
          * @brief Constructor
          */
         SetNeg();

         friend class IRegisterId<SetNeg>;

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

#endif // QUICC_ARITHMETICS_SETNEG_HPP
