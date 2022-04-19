/**
 * @file MagReynolds.hpp
 * @brief magnetic Reynolds number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_MAGREYNOLDS_HPP
#define QUICC_NONDIMENSIONAL_MAGREYNOLDS_HPP

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
    * @brief magnetic Reynolds number nondimensional number
    */
   class MagReynolds: public IRegisterId<MagReynolds>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of magnetic Reynolds number
          */
         MagReynolds(const MHDFloat value);

         friend class IRegisterId<MagReynolds>;

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

#endif // QUICC_NONDIMENSIONAL_MAGREYNOLDS_HPP
