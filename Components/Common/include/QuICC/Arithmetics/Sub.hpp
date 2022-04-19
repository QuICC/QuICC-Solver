/**
 * @file Sub.hpp
 * @brief Sub arithmetics
 */

#ifndef QUICC_ARITHMETICS_SUB_HPP
#define QUICC_ARITHMETICS_SUB_HPP

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
    * @brief Sub Arithmetics
    */
   class Sub: public IRegisterId<Sub>
   {
      public:
         /**
          * @brief Constructor
          */
         Sub();

         friend class IRegisterId<Sub>;

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

#endif // QUICC_ARITHMETICS_SUB_HPP
