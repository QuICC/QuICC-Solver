/**
 * @file I4P.hpp
 * @brief Forward transform operator Forward::I4P
 */

#ifndef QUICC_TRANSFORM_FORWARD_I4P_HPP
#define QUICC_TRANSFORM_FORWARD_I4P_HPP

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
    * @brief Forward transform operator Forward::I4P
    */
   class I4P: public IRegisterId<I4P>
   {
      public:
         /**
          * @brief Constructor
          */
         I4P();

         friend class IRegisterId<I4P>;

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

#endif // QUICC_TRANSFORM_FORWARD_I4P_HPP
