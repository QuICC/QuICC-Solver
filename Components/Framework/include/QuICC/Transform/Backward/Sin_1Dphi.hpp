/**
 * @file Sin_1Dphi.hpp
 * @brief Backward transform operator Backard::Sin_1Dphi
 */

#ifndef QUICC_TRANSFORM_BACKWARD_SIN_1DPHI_HPP
#define QUICC_TRANSFORM_BACKWARD_SIN_1DPHI_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Transform/Backward/IRegisterId.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   /**
    * @brief Backward transform operator Backard::Sin_1Dphi
    */
   class Sin_1Dphi: public IRegisterId<Sin_1Dphi>
   {
      public:
         /**
          * @brief Constructor
          */
         Sin_1Dphi();

         friend class IRegisterId<Sin_1Dphi>;

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

} // Backward
} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_BACKWARD_SIN_1DPHI_HPP
