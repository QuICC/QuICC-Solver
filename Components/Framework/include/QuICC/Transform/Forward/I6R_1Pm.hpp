/**
 * @file I6R_1Pm.hpp
 * @brief Forward transform operator Forward::I6R_1Pm
 */

#ifndef QUICC_TRANSFORM_FORWARD_I6R_1PM_HPP
#define QUICC_TRANSFORM_FORWARD_I6R_1PM_HPP

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
    * @brief Forward transform operator Forward::I6R_1Pm
    */
   class I6R_1Pm: public IRegisterId<I6R_1Pm>
   {
      public:
         /**
          * @brief Constructor
          */
         I6R_1Pm();

         friend class IRegisterId<I6R_1Pm>;

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

#endif // QUICC_TRANSFORM_FORWARD_I6R_1PM_HPP
