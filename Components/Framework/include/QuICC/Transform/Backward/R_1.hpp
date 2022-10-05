/**
 * @file R_1.hpp
 * @brief Backward transform operator Backard::R_1
 */

#ifndef QUICC_TRANSFORM_BACKWARD_R_1_HPP
#define QUICC_TRANSFORM_BACKWARD_R_1_HPP

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
    * @brief Backward transform operator Backard::R_1
    */
   class R_1: public IRegisterId<R_1>
   {
      public:
         /**
          * @brief Constructor
          */
         R_1();

         friend class IRegisterId<R_1>;

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

#endif // QUICC_TRANSFORM_BACKWARD_R_1_HPP
