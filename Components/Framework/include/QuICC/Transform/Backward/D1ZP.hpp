/**
 * @file D1ZP.hpp
 * @brief Backward transform operator Backard::D1ZP
 */

#ifndef QUICC_TRANSFORM_BACKWARD_D1ZP_HPP
#define QUICC_TRANSFORM_BACKWARD_D1ZP_HPP

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
    * @brief Backward transform operator Backard::D1ZP
    */
   class D1ZP: public IRegisterId<D1ZP>
   {
      public:
         /**
          * @brief Constructor
          */
         D1ZP();

         friend class IRegisterId<D1ZP>;

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

#endif // QUICC_TRANSFORM_BACKWARD_D1ZP_HPP
