/**
 * @file Laplh_1Sin_1.hpp
 * @brief Forward transform operator Forward::Laplh_1Sin_1
 */

#ifndef QUICC_TRANSFORM_FORWARD_LAPLH_1SIN_1_HPP
#define QUICC_TRANSFORM_FORWARD_LAPLH_1SIN_1_HPP

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
    * @brief Forward transform operator Forward::Laplh_1Sin_1
    */
   class Laplh_1Sin_1: public IRegisterId<Laplh_1Sin_1>
   {
      public:
         /**
          * @brief Constructor
          */
         Laplh_1Sin_1();

         friend class IRegisterId<Laplh_1Sin_1>;

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

#endif // QUICC_TRANSFORM_FORWARD_LAPLH_1SIN_1_HPP
