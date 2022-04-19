/**
 * @file Before.hpp
 * @brief Before SolveTiming
 */

#ifndef QUICC_SOLVETIMING_BEFORE_HPP
#define QUICC_SOLVETIMING_BEFORE_HPP

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
    * @brief Before SolveTiming
    */
   class Before: public IRegisterId<Before>
   {
      public:
         /**
          * @brief Constructor
          */
         Before();

         friend class IRegisterId<Before>;

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

#endif // QUICC_SOLVETIMING_BEFORE_HPP
