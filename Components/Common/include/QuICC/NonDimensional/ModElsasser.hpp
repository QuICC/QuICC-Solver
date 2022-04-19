/**
 * @file ModElsasser.hpp
 * @brief modified Elsasser number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_MODELSASSER_HPP
#define QUICC_NONDIMENSIONAL_MODELSASSER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/BasicTypes.hpp"
#include "QuICC/NonDimensional/IRegisterId.hpp"

namespace QuICC {

namespace NonDimensional {

   /**
    * @brief modified Elsasser number nondimensional number
    */
   class ModElsasser: public IRegisterId<ModElsasser>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of modified Elsasser number
          */
         ModElsasser(const MHDFloat value);

         friend class IRegisterId<ModElsasser>;

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

#endif // QUICC_NONDIMENSIONAL_MODELSASSER_HPP
