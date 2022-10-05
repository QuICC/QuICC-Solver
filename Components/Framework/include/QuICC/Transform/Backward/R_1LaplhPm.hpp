/**
 * @file R_1LaplhPm.hpp
 * @brief Backward transform operator Backard::R_1LaplhPm
 */

#ifndef QUICC_TRANSFORM_BACKWARD_R_1LAPLHPM_HPP
#define QUICC_TRANSFORM_BACKWARD_R_1LAPLHPM_HPP

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
    * @brief Backward transform operator Backard::R_1LaplhPm
    */
   class R_1LaplhPm: public IRegisterId<R_1LaplhPm>
   {
      public:
         /**
          * @brief Constructor
          */
         R_1LaplhPm();

         friend class IRegisterId<R_1LaplhPm>;

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

#endif // QUICC_TRANSFORM_BACKWARD_R_1LAPLHPM_HPP
