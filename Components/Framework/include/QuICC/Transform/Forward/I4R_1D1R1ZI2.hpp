/**
 * @file I4R_1D1R1ZI2.hpp
 * @brief Forward transform operator Forward::I4R_1D1R1ZI2
 */

#ifndef QUICC_TRANSFORM_FORWARD_I4R_1D1R1ZI2_HPP
#define QUICC_TRANSFORM_FORWARD_I4R_1D1R1ZI2_HPP

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
    * @brief Forward transform operator Forward::I4R_1D1R1ZI2
    */
   class I4R_1D1R1ZI2: public IRegisterId<I4R_1D1R1ZI2>
   {
      public:
         /**
          * @brief Constructor
          */
         I4R_1D1R1ZI2();

         friend class IRegisterId<I4R_1D1R1ZI2>;

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

#endif // QUICC_TRANSFORM_FORWARD_I4R_1D1R1ZI2_HPP
