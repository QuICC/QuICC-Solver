/**
 * @file DfLaplh_1.hpp
 * @brief Forward transform operator Forward::DfLaplh_1
 */

#ifndef QUICC_TRANSFORM_FORWARD_DFLAPLH_1_HPP
#define QUICC_TRANSFORM_FORWARD_DFLAPLH_1_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Transform/Forward/IRegisterId.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   /**
    * @brief Forward transform operator Forward::DfLaplh_1
    */
   class DfLaplh_1: public IRegisterId<DfLaplh_1>
   {
      public:
         /**
          * @brief Constructor
          */
         DfLaplh_1();

         friend class IRegisterId<DfLaplh_1>;

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

} // Forward
} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_FORWARD_DFLAPLH_1_HPP
