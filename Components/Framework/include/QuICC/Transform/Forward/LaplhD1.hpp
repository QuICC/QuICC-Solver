/**
 * @file LaplhD1.hpp
 * @brief Forward transform operator Forward::LaplhD1
 */

#ifndef QUICC_TRANSFORM_FORWARD_LAPLHD1_HPP
#define QUICC_TRANSFORM_FORWARD_LAPLHD1_HPP

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
    * @brief Forward transform operator Forward::LaplhD1
    */
   class LaplhD1: public IRegisterId<LaplhD1>
   {
      public:
         /**
          * @brief Constructor
          */
         LaplhD1();

         friend class IRegisterId<LaplhD1>;

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

#endif // QUICC_TRANSFORM_FORWARD_LAPLHD1_HPP
