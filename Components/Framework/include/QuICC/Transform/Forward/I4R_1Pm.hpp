/**
 * @file I4R_1Pm.hpp
 * @brief Forward transform operator Forward::I4R_1Pm
 */

#ifndef QUICC_TRANSFORM_FORWARD_I4R_1PM_HPP
#define QUICC_TRANSFORM_FORWARD_I4R_1PM_HPP

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
    * @brief Forward transform operator Forward::I4R_1Pm
    */
   class I4R_1Pm: public IRegisterId<I4R_1Pm>
   {
      public:
         /**
          * @brief Constructor
          */
         I4R_1Pm();

         friend class IRegisterId<I4R_1Pm>;

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

#endif // QUICC_TRANSFORM_FORWARD_I4R_1PM_HPP
