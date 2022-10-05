/**
 * @file R_1Pm.hpp
 * @brief Backward transform operator Backard::R_1Pm
 */

#ifndef QUICC_TRANSFORM_BACKWARD_R_1PM_HPP
#define QUICC_TRANSFORM_BACKWARD_R_1PM_HPP

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
    * @brief Backward transform operator Backard::R_1Pm
    */
   class R_1Pm: public IRegisterId<R_1Pm>
   {
      public:
         /**
          * @brief Constructor
          */
         R_1Pm();

         friend class IRegisterId<R_1Pm>;

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

#endif // QUICC_TRANSFORM_BACKWARD_R_1PM_HPP
