/**
 * @file Boundary.hpp
 * @brief Boundary ModelOperator
 */

#ifndef QUICC_MODELOPERATOR_BOUNDARY_HPP
#define QUICC_MODELOPERATOR_BOUNDARY_HPP

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
    * @brief Boundary ModelOperator
    */
   class Boundary: public IRegisterId<Boundary>
   {
      public:
         /**
          * @brief Constructor
          */
         Boundary();

         friend class IRegisterId<Boundary>;

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

#endif // QUICC_MODELOPERATOR_BOUNDARY_HPP
