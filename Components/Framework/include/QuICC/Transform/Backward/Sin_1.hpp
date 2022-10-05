/**
 * @file Sin_1.hpp
 * @brief Backward transform operator Backard::Sin_1
 */

#ifndef QUICC_TRANSFORM_BACKWARD_SIN_1_HPP
#define QUICC_TRANSFORM_BACKWARD_SIN_1_HPP

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
    * @brief Backward transform operator Backard::Sin_1
    */
   class Sin_1: public IRegisterId<Sin_1>
   {
      public:
         /**
          * @brief Constructor
          */
         Sin_1();

         friend class IRegisterId<Sin_1>;

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

#endif // QUICC_TRANSFORM_BACKWARD_SIN_1_HPP
