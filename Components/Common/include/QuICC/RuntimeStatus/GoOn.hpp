/**
 * @file GoOn.hpp
 * @brief GoOn RuntimeStatus
 */

#ifndef QUICC_RUNTIMESTATUS_GOON_HPP
#define QUICC_RUNTIMESTATUS_GOON_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/RuntimeStatus/IRegisterId.hpp"

namespace QuICC {

namespace RuntimeStatus {

   /**
    * @brief GoOn RuntimeStatus
    */
   class GoOn: public IRegisterId<GoOn>
   {
      public:
         /**
          * @brief Constructor
          */
         GoOn();

         friend class IRegisterId<GoOn>;

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

#endif // QUICC_RUNTIMESTATUS_GOON_HPP
