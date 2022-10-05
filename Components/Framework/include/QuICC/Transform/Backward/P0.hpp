/**
 * @file P0.hpp
 * @brief Backward transform operator Backard::P0
 */

#ifndef QUICC_TRANSFORM_BACKWARD_P0_HPP
#define QUICC_TRANSFORM_BACKWARD_P0_HPP

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
    * @brief Backward transform operator Backard::P0
    */
   class P0: public IRegisterId<P0>
   {
      public:
         /**
          * @brief Constructor
          */
         P0();

         friend class IRegisterId<P0>;

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

#endif // QUICC_TRANSFORM_BACKWARD_P0_HPP
