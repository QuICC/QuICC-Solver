/**
 * @file LaplhSin_1.hpp
 * @brief Forward transform operator Forward::LaplhSin_1
 */

#ifndef QUICC_TRANSFORM_FORWARD_LAPLHSIN_1_HPP
#define QUICC_TRANSFORM_FORWARD_LAPLHSIN_1_HPP

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
    * @brief Forward transform operator Forward::LaplhSin_1
    */
   class LaplhSin_1: public IRegisterId<LaplhSin_1>
   {
      public:
         /**
          * @brief Constructor
          */
         LaplhSin_1();

         friend class IRegisterId<LaplhSin_1>;

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

#endif // QUICC_TRANSFORM_FORWARD_LAPLHSIN_1_HPP
