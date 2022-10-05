/**
 * @file T.hpp
 * @brief Forward transform operator Forward::T
 */

#ifndef QUICC_TRANSFORM_FORWARD_T_HPP
#define QUICC_TRANSFORM_FORWARD_T_HPP

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
    * @brief Forward transform operator Forward::T
    */
   class T: public IRegisterId<T>
   {
      public:
         /**
          * @brief Constructor
          */
         T();

         friend class IRegisterId<T>;

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

#endif // QUICC_TRANSFORM_FORWARD_T_HPP
