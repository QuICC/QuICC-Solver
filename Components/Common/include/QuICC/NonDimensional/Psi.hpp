/**
 * @file Psi.hpp
 * @brief Psi number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_PSI_HPP
#define QUICC_NONDIMENSIONAL_PSI_HPP

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
    * @brief Psi number nondimensional number
    */
   class Psi: public IRegisterId<Psi>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Psi number
          */
         Psi(const MHDFloat value);

         friend class IRegisterId<Psi>;

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

#endif // QUICC_NONDIMENSIONAL_PSI_HPP
