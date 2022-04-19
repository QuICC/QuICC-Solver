/**
 * @file Rho.hpp
 * @brief Rho number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_RHO_HPP
#define QUICC_NONDIMENSIONAL_RHO_HPP

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
    * @brief Rho number nondimensional number
    */
   class Rho: public IRegisterId<Rho>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Rho number
          */
         Rho(const MHDFloat value);

         friend class IRegisterId<Rho>;

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

#endif // QUICC_NONDIMENSIONAL_RHO_HPP
