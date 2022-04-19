/**
 * @file ExplicitLinear.hpp
 * @brief ExplicitLinear ModelOperator
 */

#ifndef QUICC_MODELOPERATOR_EXPLICITLINEAR_HPP
#define QUICC_MODELOPERATOR_EXPLICITLINEAR_HPP

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
    * @brief ExplicitLinear ModelOperator
    */
   class ExplicitLinear: public IRegisterId<ExplicitLinear>
   {
      public:
         /**
          * @brief Constructor
          */
         ExplicitLinear();

         friend class IRegisterId<ExplicitLinear>;

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

#endif // QUICC_MODELOPERATOR_EXPLICITLINEAR_HPP
