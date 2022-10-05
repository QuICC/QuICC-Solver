/**
 * @file Laplh_1D1.hpp
 * @brief Forward transform operator Forward::Laplh_1D1
 */

#ifndef QUICC_TRANSFORM_FORWARD_LAPLH_1D1_HPP
#define QUICC_TRANSFORM_FORWARD_LAPLH_1D1_HPP

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
    * @brief Forward transform operator Forward::Laplh_1D1
    */
   class Laplh_1D1: public IRegisterId<Laplh_1D1>
   {
      public:
         /**
          * @brief Constructor
          */
         Laplh_1D1();

         friend class IRegisterId<Laplh_1D1>;

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

#endif // QUICC_TRANSFORM_FORWARD_LAPLH_1D1_HPP
