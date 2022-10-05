/**
 * @file P.hpp
 * @brief Forward transform operator Forward::P
 */

#ifndef QUICC_TRANSFORM_FORWARD_P_HPP
#define QUICC_TRANSFORM_FORWARD_P_HPP

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
    * @brief Forward transform operator Forward::P
    */
   class P: public IRegisterId<P>
   {
      public:
         /**
          * @brief Constructor
          */
         P();

         friend class IRegisterId<P>;

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

#endif // QUICC_TRANSFORM_FORWARD_P_HPP
