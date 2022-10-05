/**
 * @file Sin_1D1Sin.hpp
 * @brief Backward transform operator Backard::Sin_1D1Sin
 */

#ifndef QUICC_TRANSFORM_BACKWARD_SIN_1D1SIN_HPP
#define QUICC_TRANSFORM_BACKWARD_SIN_1D1SIN_HPP

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
    * @brief Backward transform operator Backard::Sin_1D1Sin
    */
   class Sin_1D1Sin: public IRegisterId<Sin_1D1Sin>
   {
      public:
         /**
          * @brief Constructor
          */
         Sin_1D1Sin();

         friend class IRegisterId<Sin_1D1Sin>;

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

#endif // QUICC_TRANSFORM_BACKWARD_SIN_1D1SIN_HPP
