/**
 * @file I6R_1D1R1.hpp
 * @brief Forward transform operator Forward::I6R_1D1R1
 */

#ifndef QUICC_TRANSFORM_FORWARD_I6R_1D1R1_HPP
#define QUICC_TRANSFORM_FORWARD_I6R_1D1R1_HPP

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
    * @brief Forward transform operator Forward::I6R_1D1R1
    */
   class I6R_1D1R1: public IRegisterId<I6R_1D1R1>
   {
      public:
         /**
          * @brief Constructor
          */
         I6R_1D1R1();

         friend class IRegisterId<I6R_1D1R1>;

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

#endif // QUICC_TRANSFORM_FORWARD_I6R_1D1R1_HPP
