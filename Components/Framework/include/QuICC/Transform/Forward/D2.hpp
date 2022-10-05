/**
 * @file D2.hpp
 * @brief Forward transform operator Forward::D2
 */

#ifndef QUICC_TRANSFORM_FORWARD_D2_HPP
#define QUICC_TRANSFORM_FORWARD_D2_HPP

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
    * @brief Forward transform operator Forward::D2
    */
   class D2: public IRegisterId<D2>
   {
      public:
         /**
          * @brief Constructor
          */
         D2();

         friend class IRegisterId<D2>;

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

#endif // QUICC_TRANSFORM_FORWARD_D2_HPP
