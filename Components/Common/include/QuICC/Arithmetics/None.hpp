/**
 * @file None.hpp
 * @brief None arithmetics
 */

#ifndef QUICC_ARITHMETICS_NONE_HPP
#define QUICC_ARITHMETICS_NONE_HPP

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
    * @brief None Arithmetics
    */
   class None: public IRegisterId<None>
   {
      public:
         /**
          * @brief Constructor
          */
         None();

         friend class IRegisterId<None>;

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

#endif // QUICC_ARITHMETICS_NONE_HPP
