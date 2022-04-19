/**
 * @file After.hpp
 * @brief After SolveTiming
 */

#ifndef QUICC_SOLVETIMING_AFTER_HPP
#define QUICC_SOLVETIMING_AFTER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/SolveTiming/IRegisterId.hpp"

namespace QuICC {

namespace SolveTiming {

   /**
    * @brief After SolveTiming
    */
   class After: public IRegisterId<After>
   {
      public:
         /**
          * @brief Constructor
          */
         After();

         friend class IRegisterId<After>;

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

#endif // QUICC_SOLVETIMING_AFTER_HPP
