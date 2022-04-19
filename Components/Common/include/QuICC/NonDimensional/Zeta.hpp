/**
 * @file Zeta.hpp
 * @brief Zeta number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_ZETA_HPP
#define QUICC_NONDIMENSIONAL_ZETA_HPP

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
    * @brief Zeta number nondimensional number
    */
   class Zeta: public IRegisterId<Zeta>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Zeta number
          */
         Zeta(const MHDFloat value);

         friend class IRegisterId<Zeta>;

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

#endif // QUICC_NONDIMENSIONAL_ZETA_HPP
