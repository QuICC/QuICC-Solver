/**
 * @file Stop.hpp
 * @brief Stop RuntimeStatus
 */

#ifndef QUICC_RUNTIMESTATUS_STOP_HPP
#define QUICC_RUNTIMESTATUS_STOP_HPP

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
    * @brief Stop RuntimeStatus
    */
   class Stop: public IRegisterId<Stop>
   {
      public:
         /**
          * @brief Constructor
          */
         Stop();

         friend class IRegisterId<Stop>;

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

#endif // QUICC_RUNTIMESTATUS_STOP_HPP
