/**
 * @file LaplhSin_1Dphi.hpp
 * @brief Forward transform operator Forward::LaplhSin_1Dphi
 */

#ifndef QUICC_TRANSFORM_FORWARD_LAPLHSIN_1DPHI_HPP
#define QUICC_TRANSFORM_FORWARD_LAPLHSIN_1DPHI_HPP

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
    * @brief Forward transform operator Forward::LaplhSin_1Dphi
    */
   class LaplhSin_1Dphi: public IRegisterId<LaplhSin_1Dphi>
   {
      public:
         /**
          * @brief Constructor
          */
         LaplhSin_1Dphi();

         friend class IRegisterId<LaplhSin_1Dphi>;

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

#endif // QUICC_TRANSFORM_FORWARD_LAPLHSIN_1DPHI_HPP
