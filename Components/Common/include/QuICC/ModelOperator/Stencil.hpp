/**
 * @file Stencil.hpp
 * @brief Stencil ModelOperator
 */

#ifndef QUICC_MODELOPERATOR_STENCIL_HPP
#define QUICC_MODELOPERATOR_STENCIL_HPP

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
    * @brief Stencil ModelOperator
    */
   class Stencil: public IRegisterId<Stencil>
   {
      public:
         /**
          * @brief Constructor
          */
         Stencil();

         friend class IRegisterId<Stencil>;

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

#endif // QUICC_MODELOPERATOR_STENCIL_HPP
