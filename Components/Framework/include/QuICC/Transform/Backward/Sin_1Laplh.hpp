/**
 * @file Sin_1Laplh.hpp
 * @brief Backward transform operator Backard::Sin_1Laplh
 */

#ifndef QUICC_TRANSFORM_BACKWARD_SIN_1LAPLH_HPP
#define QUICC_TRANSFORM_BACKWARD_SIN_1LAPLH_HPP

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
    * @brief Backward transform operator Backard::Sin_1Laplh
    */
   class Sin_1Laplh: public IRegisterId<Sin_1Laplh>
   {
      public:
         /**
          * @brief Constructor
          */
         Sin_1Laplh();

         friend class IRegisterId<Sin_1Laplh>;

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

#endif // QUICC_TRANSFORM_BACKWARD_SIN_1LAPLH_HPP
