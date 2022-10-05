/**
 * @file Laplh2.hpp
 * @brief Forward transform operator Forward::Laplh2
 */

#ifndef QUICC_TRANSFORM_FORWARD_LAPLH2_HPP
#define QUICC_TRANSFORM_FORWARD_LAPLH2_HPP

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
    * @brief Forward transform operator Forward::Laplh2
    */
   class Laplh2: public IRegisterId<Laplh2>
   {
      public:
         /**
          * @brief Constructor
          */
         Laplh2();

         friend class IRegisterId<Laplh2>;

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

#endif // QUICC_TRANSFORM_FORWARD_LAPLH2_HPP
