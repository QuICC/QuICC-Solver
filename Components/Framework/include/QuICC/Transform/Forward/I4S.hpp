/**
 * @file I4S.hpp
 * @brief Forward transform operator Forward::I4S
 */

#ifndef QUICC_TRANSFORM_FORWARD_I4S_HPP
#define QUICC_TRANSFORM_FORWARD_I4S_HPP

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
    * @brief Forward transform operator Forward::I4S
    */
   class I4S: public IRegisterId<I4S>
   {
      public:
         /**
          * @brief Constructor
          */
         I4S();

         friend class IRegisterId<I4S>;

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

#endif // QUICC_TRANSFORM_FORWARD_I4S_HPP
