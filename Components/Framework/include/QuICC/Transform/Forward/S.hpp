/**
 * @file S.hpp
 * @brief Forward transform operator Forward::S
 */

#ifndef QUICC_TRANSFORM_FORWARD_S_HPP
#define QUICC_TRANSFORM_FORWARD_S_HPP

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
    * @brief Forward transform operator Forward::S
    */
   class S: public IRegisterId<S>
   {
      public:
         /**
          * @brief Constructor
          */
         S();

         friend class IRegisterId<S>;

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

#endif // QUICC_TRANSFORM_FORWARD_S_HPP
