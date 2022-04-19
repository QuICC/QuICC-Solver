/**
 * @file Kappa.hpp
 * @brief Kappa number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_KAPPA_HPP
#define QUICC_NONDIMENSIONAL_KAPPA_HPP

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
    * @brief Kappa number nondimensional number
    */
   class Kappa: public IRegisterId<Kappa>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Kappa number
          */
         Kappa(const MHDFloat value);

         friend class IRegisterId<Kappa>;

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

#endif // QUICC_NONDIMENSIONAL_KAPPA_HPP
