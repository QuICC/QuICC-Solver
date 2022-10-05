/**
 * @file D1.hpp
 * @brief Forward transform operator Forward::D1
 */

#ifndef QUICC_TRANSFORM_FORWARD_D1_HPP
#define QUICC_TRANSFORM_FORWARD_D1_HPP

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
    * @brief Forward transform operator Forward::D1
    */
   class D1: public IRegisterId<D1>
   {
      public:
         /**
          * @brief Constructor
          */
         D1();

         friend class IRegisterId<D1>;

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

#endif // QUICC_TRANSFORM_FORWARD_D1_HPP
