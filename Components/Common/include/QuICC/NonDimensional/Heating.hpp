/**
 * @file Heating.hpp
 * @brief Heating number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_HEATING_HPP
#define QUICC_NONDIMENSIONAL_HEATING_HPP

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
    * @brief Heating number nondimensional number
    */
   class Heating: public IRegisterId<Heating>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Heating number
          */
         Heating(const MHDFloat value);

         friend class IRegisterId<Heating>;

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

#endif // QUICC_NONDIMENSIONAL_HEATING_HPP
