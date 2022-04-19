/**
 * @file ExplicitNonlinear.hpp
 * @brief ExplicitNonlinear ModelOperator
 */

#ifndef QUICC_MODELOPERATOR_EXPLICITNONLINEAR_HPP
#define QUICC_MODELOPERATOR_EXPLICITNONLINEAR_HPP

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
    * @brief ExplicitNonlinear ModelOperator
    */
   class ExplicitNonlinear: public IRegisterId<ExplicitNonlinear>
   {
      public:
         /**
          * @brief Constructor
          */
         ExplicitNonlinear();

         friend class IRegisterId<ExplicitNonlinear>;

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

#endif // QUICC_MODELOPERATOR_EXPLICITNONLINEAR_HPP
