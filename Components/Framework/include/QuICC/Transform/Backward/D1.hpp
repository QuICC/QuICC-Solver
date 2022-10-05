/**
 * @file D1.hpp
 * @brief Backward transform operator Backward::D1
 */

#ifndef QUICC_TRANSFORM_BACKWARD_D1_HPP
#define QUICC_TRANSFORM_BACKWARD_D1_HPP

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
    * @brief Backward transform operator Backward::D1
    */
   class D1: public IRegisterId<D1>
   {
      public:
         /**
          * @brief Constructor
          */
         D1();

         friend class IRegisterId<D1>;

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

#endif // QUICC_TRANSFORM_BACKWARD_D1_HPP
