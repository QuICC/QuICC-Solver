/**
 * @file D3.hpp
 * @brief Backward transform operator Backard::D3
 */

#ifndef QUICC_TRANSFORM_BACKWARD_D3_HPP
#define QUICC_TRANSFORM_BACKWARD_D3_HPP

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
    * @brief Backward transform operator Backard::D3
    */
   class D3: public IRegisterId<D3>
   {
      public:
         /**
          * @brief Constructor
          */
         D3();

         friend class IRegisterId<D3>;

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

#endif // QUICC_TRANSFORM_BACKWARD_D3_HPP
