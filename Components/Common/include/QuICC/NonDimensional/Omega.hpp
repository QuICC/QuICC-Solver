/**
 * @file Omega.hpp
 * @brief Omega number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_OMEGA_HPP
#define QUICC_NONDIMENSIONAL_OMEGA_HPP

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
    * @brief Omega number nondimensional number
    */
   class Omega: public IRegisterId<Omega>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Omega number
          */
         Omega(const MHDFloat value);

         friend class IRegisterId<Omega>;

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

#endif // QUICC_NONDIMENSIONAL_OMEGA_HPP
