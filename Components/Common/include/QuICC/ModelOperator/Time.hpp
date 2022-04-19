/**
 * @file Time.hpp
 * @brief Time ModelOperator
 */

#ifndef QUICC_MODELOPERATOR_TIME_HPP
#define QUICC_MODELOPERATOR_TIME_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/ModelOperator/IRegisterId.hpp"

namespace QuICC {

namespace ModelOperator {

   /**
    * @brief Time ModelOperator
    */
   class Time: public IRegisterId<Time>
   {
      public:
         /**
          * @brief Constructor
          */
         Time();

         friend class IRegisterId<Time>;

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

#endif // QUICC_MODELOPERATOR_TIME_HPP
