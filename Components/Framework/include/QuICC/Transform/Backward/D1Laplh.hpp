/**
 * @file D1Laplh.hpp
 * @brief Backward transform operator Backard::D1Laplh
 */

#ifndef QUICC_TRANSFORM_BACKWARD_D1LAPLH_HPP
#define QUICC_TRANSFORM_BACKWARD_D1LAPLH_HPP

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
    * @brief Backward transform operator Backard::D1Laplh
    */
   class D1Laplh: public IRegisterId<D1Laplh>
   {
      public:
         /**
          * @brief Constructor
          */
         D1Laplh();

         friend class IRegisterId<D1Laplh>;

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

#endif // QUICC_TRANSFORM_BACKWARD_D1LAPLH_HPP
