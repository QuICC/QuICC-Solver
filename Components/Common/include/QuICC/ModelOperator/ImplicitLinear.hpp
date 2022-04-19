/**
 * @file ImplicitLinear.hpp
 * @brief ImplicitLinear ModelOperator
 */

#ifndef QUICC_MODELOPERATOR_IMPLICITLINEAR_HPP
#define QUICC_MODELOPERATOR_IMPLICITLINEAR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/ModelOperator/IRegisterId.hpp"

namespace QuICC {

namespace ModelOperator {

   /**
    * @brief ImplicitLinear ModelOperator
    */
   class ImplicitLinear: public IRegisterId<ImplicitLinear>
   {
      public:
         /**
          * @brief Constructor
          */
         ImplicitLinear();

         friend class IRegisterId<ImplicitLinear>;

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

}
}

#endif // QUICC_MODELOPERATOR_IMPLICITLINEAR_HPP
