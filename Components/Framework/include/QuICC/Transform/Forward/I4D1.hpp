/**
 * @file I4D1.hpp
 * @brief Forward transform operator Forward::I4D1
 */

#ifndef QUICC_TRANSFORM_FORWARD_I4D1_HPP
#define QUICC_TRANSFORM_FORWARD_I4D1_HPP

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
    * @brief Forward transform operator Forward::I4D1
    */
   class I4D1: public IRegisterId<I4D1>
   {
      public:
         /**
          * @brief Constructor
          */
         I4D1();

         friend class IRegisterId<I4D1>;

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

#endif // QUICC_TRANSFORM_FORWARD_I4D1_HPP
