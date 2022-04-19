/**
 * @file ExplicitNextstep.hpp
 * @brief ExplicitNextstep ModelOperator
 */

#ifndef QUICC_MODELOPERATOR_EXPLICITNEXTSTEP_HPP
#define QUICC_MODELOPERATOR_EXPLICITNEXTSTEP_HPP

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
    * @brief ExplicitNextstep ModelOperator
    */
   class ExplicitNextstep: public IRegisterId<ExplicitNextstep>
   {
      public:
         /**
          * @brief Constructor
          */
         ExplicitNextstep();

         friend class IRegisterId<ExplicitNextstep>;

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

#endif // QUICC_MODELOPERATOR_EXPLICITNEXTSTEP_HPP
