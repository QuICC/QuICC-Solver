/**
 * @file R_2.hpp
 * @brief Backward transform operator Backard::R_2
 */

#ifndef QUICC_TRANSFORM_BACKWARD_R_2_HPP
#define QUICC_TRANSFORM_BACKWARD_R_2_HPP

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
    * @brief Backward transform operator Backard::R_2
    */
   class R_2: public IRegisterId<R_2>
   {
      public:
         /**
          * @brief Constructor
          */
         R_2();

         friend class IRegisterId<R_2>;

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

#endif // QUICC_TRANSFORM_BACKWARD_R_2_HPP
