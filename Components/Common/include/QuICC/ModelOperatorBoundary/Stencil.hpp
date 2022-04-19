/**
 * @file Stencil.hpp
 * @brief Stencil ModelOperatorBoundary
 */

#ifndef QUICC_MODELOPERATORBOUNDARY_STENCIL_HPP
#define QUICC_MODELOPERATORBOUNDARY_STENCIL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/ModelOperatorBoundary/IRegisterId.hpp"

namespace QuICC {

namespace ModelOperatorBoundary {

   /**
    * @brief Stencil ModelOperatorBoundary
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

#endif // QUICC_MODELOPERATORBOUNDARY_STENCIL_HPP
