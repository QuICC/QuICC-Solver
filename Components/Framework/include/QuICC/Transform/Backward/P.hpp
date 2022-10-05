/**
 * @file P.hpp
 * @brief Backward transform operator Backard::P
 */

#ifndef QUICC_TRANSFORM_BACKWARD_P_HPP
#define QUICC_TRANSFORM_BACKWARD_P_HPP

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
    * @brief Backward transform operator Backard::P
    */
   class P: public IRegisterId<P>
   {
      public:
         /**
          * @brief Constructor
          */
         P();

         friend class IRegisterId<P>;

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

#endif // QUICC_TRANSFORM_BACKWARD_P_HPP
