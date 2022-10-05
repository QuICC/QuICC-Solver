/**
 * @file R_1D1R1.hpp
 * @brief Backward transform operator Backard::R_1D1R1
 */

#ifndef QUICC_TRANSFORM_BACKWARD_R_1D1R1_HPP
#define QUICC_TRANSFORM_BACKWARD_R_1D1R1_HPP

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
    * @brief Backward transform operator Backard::R_1D1R1
    */
   class R_1D1R1: public IRegisterId<R_1D1R1>
   {
      public:
         /**
          * @brief Constructor
          */
         R_1D1R1();

         friend class IRegisterId<R_1D1R1>;

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

#endif // QUICC_TRANSFORM_BACKWARD_R_1D1R1_HPP
